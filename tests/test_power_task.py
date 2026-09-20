"""Run the production task loop with a deterministic clock on the host."""

import pathlib
import re
import subprocess
import tempfile
import unittest


class PowerTaskTest(unittest.TestCase):
    def test_periods_order_delay_wraparound_and_no_supercap(self):
        module = pathlib.Path(__file__).resolve().parents[1]
        source = (module / "PowerControl.hpp").read_text()
        signature = "static void ThreadFunc(PowerControl* self)"
        self.assertIn(signature, source, "PowerControl must own the power task")
        start = source.index(signature)
        body = source.index("{", start)
        depth = 1
        end = body + 1
        while depth:
            depth += (source[end] == "{") - (source[end] == "}")
            end += 1
        loop = source[start:end]
        super_source = (module.parent / "SuperPower/SuperPower.hpp").read_text()
        constants = "\n".join(re.findall(
            r"static constexpr uint32_t (?:BACKGROUND|COMMAND)_PERIOD_MS = .*?;",
            source + super_source,
        ))
        harness = r'''
#include <cassert>
#include <cstdint>
#include <utility>
#include <vector>
static std::vector<uint32_t> times;
static std::vector<std::pair<char, uint32_t>> events;
static unsigned index_now;
struct Done {};
namespace LibXR {
struct Timebase {
  static uint32_t GetMilliseconds() { return times.at(index_now); }
};
struct Thread {
  static void SleepUntil(uint32_t&, uint32_t period) {
    assert(period == 1);
    if (++index_now == times.size()) throw Done{};
  }
};
}
struct SuperPower {
  CONSTANTS
  void Update() { events.emplace_back('S', times.at(index_now)); }
};
struct PowerControl {
  CONSTANTS
  SuperPower* superpower_;
  LibXR::Thread thread_;
  void BackgroundUpdate() { events.emplace_back('P', times.at(index_now)); }
  LOOP
};
static void Run(std::vector<uint32_t> clock, bool cap = true) {
  times = std::move(clock);
  index_now = 0;
  events.clear();
  SuperPower super;
  PowerControl power{cap ? &super : nullptr, {}};
  try { PowerControl::ThreadFunc(&power); } catch (const Done&) {}
}
int main() {
  Run({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
  const std::vector<std::pair<char, uint32_t>> NORMAL{
    {'S',0}, {'P',0}, {'P',1}, {'P',2}, {'P',3}, {'P',4},
    {'S',5}, {'P',5}, {'P',6}, {'P',7}, {'P',8}, {'P',9},
    {'S',10}, {'P',10}};
  assert(events == NORMAL);
  Run({0, 1, 12, 12, 13, 17});
  const std::vector<std::pair<char, uint32_t>> DELAYED{
    {'S',0}, {'P',0}, {'P',1}, {'S',12}, {'P',12},
    {'P',12}, {'P',13}, {'S',17}, {'P',17}};
  assert(events == DELAYED);
  Run({UINT32_MAX - 2, UINT32_MAX, 0, 1, 2});
  const std::vector<std::pair<char, uint32_t>> WRAPPED{
    {'S',UINT32_MAX - 2}, {'P',UINT32_MAX - 2}, {'P',UINT32_MAX},
    {'P',0}, {'P',1}, {'S',2}, {'P',2}};
  assert(events == WRAPPED);
  Run({0, 1, 5}, false);
  const std::vector<std::pair<char, uint32_t>> NO_CAP{
    {'P',0}, {'P',1}, {'P',5}};
  assert(events == NO_CAP);
}
'''.replace("CONSTANTS", constants).replace("LOOP", loop)
        with tempfile.TemporaryDirectory() as directory:
            cpp = pathlib.Path(directory) / "test.cpp"
            binary = pathlib.Path(directory) / "test"
            cpp.write_text(harness)
            subprocess.run(["c++", "-std=c++20", "-Wall", "-Wextra", "-Werror",
                            str(cpp), "-o", str(binary)], check=True)
            subprocess.run([str(binary)], check=True)


if __name__ == "__main__":
    unittest.main()
