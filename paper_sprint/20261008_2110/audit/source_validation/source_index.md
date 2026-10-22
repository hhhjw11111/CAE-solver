# 精确源码索引

路径以 solver_v2 为根；真实文件绝对路径同时写入 source_index.json。行号为本轮文件中的1基行号。

| 路径与行 | 内容 | SHA-256 |
|---|---|---|
| `src/core/field_layout_m.f90:5–15` | 七场排序 | `0b02e79e04674d322966a1076e64dbb944ffb1634d6d007a7e24e7acde335c96` |
| `src/core/constants_m.f90:6–7` | q与kB | `dc6872db24efbbb10627ef53a0b46fb182a9b61940a66c2497456d99c31c0d15` |
| `src/physics/zno_candidate_m.f90:8–45` | 全部ZnO候选常数与模型开关 | `407519d9ec931d9d29ed1a32d2204f34a74238f1ed0ccb2cbeaf3b2c63f81026` |
| `src/physics/transport_m.f90:40–58` | N0冻结及NT温导数 | `7215021b45f41f4e0cb7654d1bbe1318275581af2eaade621b9d0cd5dee48bd1` |
| `src/physics/thermo_piezo_m.f90:11–24` | 热压电本构 | `0acecab76b6e6468b90107c75d32b69b74a08fc0599e61584a3cf5f71ea0e65d` |
| `src/physics/material_point_m.f90:107–149` | 密度、载流子电流、热源 | `649601bea7dd70704f9eecfe3d4b28de7cc80cb7b33c4ab6a8ddb8568e2ee0aa` |
| `src/physics/material_point_m.f90:151–204` | 一致材料导数及LQ替换 | `649601bea7dd70704f9eecfe3d4b28de7cc80cb7b33c4ab6a8ddb8568e2ee0aa` |
| `src/physics/recombination_m.f90:10–18` | R=0 | `06133b7b899ca05634e56b40afc9fe3f7509d3b15033a340b30d9a88200cf800` |
| `src/fem/element_kernel_m.f90:63–140` | 积分、生产残差/切线和格式入口 | `ce69b2a924e58a91f2edde78bd0568fc3689371ad9a3e65dd4c30f70db0fbda8` |
| `src/fem/element_kernel_m.f90:143–176` | log-weighted测试及其状态导数 | `ce69b2a924e58a91f2edde78bd0568fc3689371ad9a3e65dd4c30f70db0fbda8` |
| `src/fem/element_kernel_m.f90:178–236` | 普通弱式状态和测试算子 | `ce69b2a924e58a91f2edde78bd0568fc3689371ad9a3e65dd4c30f70db0fbda8` |
| `src/app/f1_case_runner_m.f90:157–178` | LQ三阶段 | `47fe98f1b659509ffcf0cfd1aec8141c7d5ea44aa61617a56b3bee0be010a6c6` |
| `src/app/f1_case_runner_m.f90:315–329` | 原验收门槛 | `47fe98f1b659509ffcf0cfd1aec8141c7d5ea44aa61617a56b3bee0be010a6c6` |
| `src/app/f1_case_runner_m.f90:543–595` | 全域电极与A/B/C | `47fe98f1b659509ffcf0cfd1aec8141c7d5ea44aa61617a56b3bee0be010a6c6` |
| `src/app/f1_case_runner_m.f90:598–657` | 半域对称与刚体约束 | `47fe98f1b659509ffcf0cfd1aec8141c7d5ea44aa61617a56b3bee0be010a6c6` |
| `src/app/f1_case_runner_m.f90:659–680` | LQ固定字段 | `47fe98f1b659509ffcf0cfd1aec8141c7d5ea44aa61617a56b3bee0be010a6c6` |
| `src/io/native_input_m.f90:59–62` | 新格式默认关闭 | `ed1b5979bf08fb8a9dd4261a8868e9dd40bbcbed87c3bd68b27ad3e2e72b8451` |
| `src/verification/verification_driver_m.f90:41–76` | 真实verify入口 | `ab204284f55a2924ca8e48f9520cf3e75bae2a7d1ef45c025048070b7b1c8de9` |
| `src/app/main.f90:29–90` | version/verify/run-f1命令分发 | `acec6d3f0674d523205d22b5d0dfeb2cd95b6416d9c8e37945b9d78a52185fc8` |
| `tests/hv_log_weighted_production_bar.f90:89–113` | 生产核条带采样误差定义 | `fbb1472fd7ce1b416bb208bf24ae34f660c8d564acc1861cb16c9c9234a3d189` |
| `tests/hv_log_weighted_production_bar.f90:123–153` | 自定义缩约装配调用生产单元 | `fbb1472fd7ce1b416bb208bf24ae34f660c8d564acc1861cb16c9c9234a3d189` |
| `tests/hv_log_weighted_production_bar.f90:187–213` | 实际Pe/网格及容差 | `fbb1472fd7ce1b416bb208bf24ae34f660c8d564acc1861cb16c9c9234a3d189` |
