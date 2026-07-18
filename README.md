# PowerControl

`PowerControl` 延续原模块的两头文件结构：

- `RLS.hpp`：固定维度、固定存储的递归最小二乘估计器。
- `PowerControl.hpp`：功率模型、能量预算、共享分配和模块接口。

没有独立算法层、后台线程或堆分配。底盘控制周期直接提交反馈和请求，
再同步调用一次 `OutputLimit()`。

## 控制流程

1. `SetMotorFeedback3508/6020` 提交实际电流指令、实际转速和可选在线掩码，
   供下一条新的超电功率遥测进行 RLS 更新。
2. `SetMotorData3508/6020` 提交本周期请求电流、转速、跟踪误差和同一在线
   掩码；离线电机不参与分配且输出为零。
3. `SetPowerRequest` 或 `SetBoostRequested` 选择 `Normal` / `Boost`。
4. `OutputLimit()` 读取一次遥测，更新模型与能量预算，并对全部 M3508、
   GM6020 执行一次共享功率分配。
5. `GetPowerControlData()` 返回限幅电流和诊断状态。

功率模型保持 RM2024 形式：

```text
P = tau * omega + k1 * abs(omega) + k2 * tau^2
```

RLS 只消费反馈数据，每个 `chassis_power_sequence` 最多尝试一次；任一配置
电机离线、激励不足或单次创新异常时整次样本都会被拒绝。在 1000 ms
有效期内累计 32 个合格的正向大创新样本后，模块以当前可信参数重置协方差并进入
受限恢复更新；短暂超电离线或无效样本不会立即清空证据，超过有效期才重置。
GM6020 的固定模型功率和底盘静态损耗会先从实测功率中扣除。

## 内部预算

模块内部固定使用 RM2024 的基础/满能量双环，并直接持有两个
`LibXR::PID<float>`。两者参数均为 `P=50`、`I=0`、`D=0.2`；
`LibXR::Timebase` 只提供控制周期 `dt`，首个样本使用零外部导数，避免启动
微分尖峰；超电原始能量与裁判缓冲能量之间切换时也用零导数重新播种，避免
跨单位导数。PID 输出限制为超电额外功率范围 `300 W`。`Normal` 请求裁判
功率上限，`Boost` 请求由当前能量状态允许的上限。超电不可用时改用裁判缓冲
能量；两者都不可用时复位两个控制器并进入保守离线降级。来源恢复后按当前
遥测立即计算预算（immediate recovery），不存在按调用次数增长的恢复状态。

冷启动尚无有效裁判上限时固定采用最低支持上限 `45 W`，不会高于一级步兵
功率等级。

裁判 `0x0201` 功率上限和 `0x0202` 缓冲能量的新鲜度由 `SuperPower` 集中
维护，任一兼容总标志都不会把已过期的单项重新判为在线。`0x0201` 新鲜度
过期后仍沿用最后一个有效上限；这是明确的控制策略。

## 配置

```yaml
- id: power_control
  name: PowerControl
  constructor_args:
    superpower: '@&superpower'
    chassis_static_power_loss: 4.5
    motor_count_3508: 4
    motor_count_6020: 0
```

模型系数、RLS 边界、能量目标、离线比例、功率余量和 PD 参数均为
`PowerControl.hpp` 内部安全默认值，正常接入无需展开调参表。
