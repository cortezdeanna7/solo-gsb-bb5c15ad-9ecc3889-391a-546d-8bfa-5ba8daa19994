# DNA Restriction Digest Simulator

基于 Python + FastAPI 的 DNA 限制性内切酶平切模拟服务。接收 ACGT 序列、线性/环状标记和酶规则，搜索识别串及其反向互补串（保留重叠匹配），合并所有切点同时切割，返回各酶切点与酶切片段。亦可按每种酶的切开概率做部分酶切，精确计算各长度片段的期望数量；还可为多个 DNA 样本设计最小酶组合（电泳面板），使任意两个样本的条带图谱都可区分。

## 启动

单条命令（需要 Docker 与 Docker Compose）：

```bash
docker compose up --build
```

镜像构建时安装依赖，容器启动时会再次执行 `pip install -r requirements.txt` 确认依赖就绪，然后启动 uvicorn。

## 配置

| 环境变量 | 默认值 | 说明 |
|---|---|---|
| `PORT` | `8000` | 宿主机映射端口，如 `PORT=9000 docker compose up --build` |

## API 地址

- 服务地址： `http://localhost:8000`
- 交互式文档（Swagger UI）： `http://localhost:8000/docs`
- 健康检查： `GET http://localhost:8000/health` → `{"status": "ok"}`
- 酶切模拟： `POST http://localhost:8000/digest`
- 部分酶切（概率）期望： `POST http://localhost:8000/partial-digest`
- 酶面板设计： `POST http://localhost:8000/design-panel`

## 接口说明

### `POST /digest`

请求体：

| 字段 | 类型 | 说明 |
|---|---|---|
| `sequence` | string | DNA 序列，仅含 A/C/G/T（大小写均可，内部统一大写） |
| `topology` | string | `"linear"`（默认）或 `"circular"` |
| `enzymes` | array | 酶规则列表，每条含 `name`（唯一）、`site`（IUPAC 识别串）、`cut_offset`（串内切割偏移 k，0 ≤ k ≤ 串长） |

规则要点：

- 同时搜索识别串及其反向互补串，保留重叠匹配；正向切点 = 起点 + k，反向切点 = 起点 + 串长 − k。
- 环状序列允许跨首尾的匹配，切点按序列长度取模；识别串不得长于序列（否则返回 422）。
- 切点为零基坐标，同酶同点去重；所有酶的切点合并后同时切割。
- 片段按起点排序，含 `start`、`length` 和原序列方向的碱基串；环状末段绕回首段。
- 无切点时返回原序列作为唯一片段；不输出空片段。

示例（线性，双酶切）：

```bash
curl -s -X POST http://localhost:8000/digest \
  -H 'Content-Type: application/json' \
  -d '{
        "sequence": "AGGATCCGAATTCA",
        "topology": "linear",
        "enzymes": [
          {"name": "BamHI", "site": "GGATCC", "cut_offset": 1},
          {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}
        ]
      }'
```

响应：

```json
{
  "cuts": [
    {"name": "BamHI", "cuts": [2, 6]},
    {"name": "EcoRI", "cuts": [8, 12]}
  ],
  "fragments": [
    {"start": 0,  "length": 2, "sequence": "AG"},
    {"start": 2,  "length": 4, "sequence": "GATC"},
    {"start": 6,  "length": 2, "sequence": "CG"},
    {"start": 8,  "length": 4, "sequence": "AATT"},
    {"start": 12, "length": 2, "sequence": "CA"}
  ]
}
```

环状示例（末段绕回首段）：

```bash
curl -s -X POST http://localhost:8000/digest \
  -H 'Content-Type: application/json' \
  -d '{
        "sequence": "GAATTC",
        "topology": "circular",
        "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}]
      }'
```

```json
{
  "cuts": [{"name": "EcoRI", "cuts": [1, 5]}],
  "fragments": [
    {"start": 1, "length": 4, "sequence": "AATT"},
    {"start": 5, "length": 2, "sequence": "CG"}
  ]
}
```

### `POST /partial-digest`

部分酶切：每种酶在每个切点以概率 `p` 独立切开，精确计算片段长度的期望分布（不采样、不枚举切口子集，以有效切点数 K 计时间复杂度 O(K²)）。

请求体：

| 字段 | 类型 | 说明 |
|---|---|---|
| `sequence` | string | DNA 序列，仅含 A/C/G/T（大小写均可），最长 2000 bp，超长返回 422 |
| `topology` | string | `"linear"`（默认）或 `"circular"` |
| `enzymes` | array | 酶规则列表，每条含 `name`（唯一）、`site`（IUPAC 识别串）、`cut_offset`（0 ≤ k ≤ 串长）及 `p`（切开概率，0～1，非法返回 422）；空列表视为不切 |

规则要点：

- 匹配与切点去重规则与 `/digest` 完全一致（识别串 + 反向互补、保留重叠、同酶同点去重）。
- 每种酶在各切点独立以概率 `p` 切开；同一位置被多种酶切中时，任一切开即形成切口，该位置切口概率为 `1 − ∏(1 − p_i)`；不同位置的切口事件相互独立。
- 线性分子两端（位置 0 与序列末尾）不算随机切口，不产生片段也不影响无切口概率。
- 响应 `fragments` 按长度升序给出片段期望数量：同长度累加、可为小数、期望为零的项省略。
- `no_cut_probability` 为整条分子无任何切口的概率。环状分子无切口或恰有一个切口时都贡献一个全长片段，但 `no_cut_probability` 只计无切口这一种情况。
- `p` 支持 0 与 1：`p = 0` 的酶等同不切，`p = 1` 时退化为确定性酶切结果。

示例（线性，EcoRI 每点 50% 概率切开）：

```bash
curl -s -X POST http://localhost:8000/partial-digest \
  -H 'Content-Type: application/json' \
  -d '{
        "sequence": "ACGTGAATTCGGGAATTCAA",
        "topology": "linear",
        "enzymes": [
          {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1, "p": 0.5}
        ]
      }'
```

响应：

```json
{
  "fragments": [
    {"length": 3,  "expected_count": 0.5},
    {"length": 4,  "expected_count": 0.75},
    {"length": 5,  "expected_count": 0.5},
    {"length": 7,  "expected_count": 0.25},
    {"length": 8,  "expected_count": 0.25},
    {"length": 9,  "expected_count": 0.25},
    {"length": 11, "expected_count": 0.125},
    {"length": 12, "expected_count": 0.0625},
    {"length": 13, "expected_count": 0.125},
    {"length": 15, "expected_count": 0.0625},
    {"length": 17, "expected_count": 0.0625},
    {"length": 20, "expected_count": 0.0625}
  ],
  "no_cut_probability": 0.0625
}
```

### `POST /design-panel`

为 2～12 个 DNA 样本挑选能区分所有样本对的最少酶组合（电泳面板）。

请求体：

| 字段 | 类型 | 说明 |
|---|---|---|
| `samples` | array | 2～12 个样本，每个含 `name`（唯一）、`sequence`（仅 A/C/G/T，大小写均可）、`topology`（`"linear"` 默认 / `"circular"`） |
| `enzymes` | array | 1～12 种候选酶，格式与 `/digest` 相同：`name`（唯一）、`site`（IUPAC）、`cut_offset`（0 ≤ k ≤ 串长，且识别串不得长于任一样本序列） |
| `tolerance_bp` | int | 条带容差（bp），非负整数 |

规则要点：

- 复用 `/digest` 的酶切逻辑：每个样本分别用每种候选酶做单酶切，每种酶独立成一条泳道。
- 泳道条带 = 片段长度升序列表，重复条带保留。
- 两条泳道条带数量相同且对应长度差均 ≤ `tolerance_bp` 时视为不可区分，否则该酶可区分这对样本。
- 选择能区分全部样本对的最少酶组；同数量方案按排序后的酶名列表字典序取最小。
- 无解（全部候选酶合起来仍无法区分某些样本对）时：`enzymes`、`lanes`、`pairs` 为空，`undistinguishable_pairs` 列出这些样本对。

响应字段：

| 字段 | 说明 |
|---|---|
| `enzymes` | 所选酶名（排序后）；无解时为空 |
| `lanes` | 每个样本在所选泳道中的条带（`bands`：酶名 → 升序长度列表），按样本输入顺序 |
| `pairs` | 每对样本及可区分它们的所选酶 |
| `undistinguishable_pairs` | 全部候选酶仍无法区分的样本对（有解时为空） |

示例：

```bash
curl -s -X POST http://localhost:8000/design-panel \
  -H 'Content-Type: application/json' \
  -d '{
        "samples": [
          {"name": "wild-type", "sequence": "AAGAATTCTT"},
          {"name": "mutant",    "sequence": "AAGGATCCTT"}
        ],
        "enzymes": [
          {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1},
          {"name": "BamHI", "site": "GGATCC", "cut_offset": 1}
        ],
        "tolerance_bp": 0
      }'
```

响应（两种酶单独都能区分该样本对，按字典序取 BamHI）：

```json
{
  "enzymes": ["BamHI"],
  "lanes": [
    {"sample": "wild-type", "bands": {"BamHI": [10]}},
    {"sample": "mutant",    "bands": {"BamHI": [3, 3, 4]}}
  ],
  "pairs": [
    {"samples": ["wild-type", "mutant"], "enzymes": ["BamHI"]}
  ],
  "undistinguishable_pairs": []
}
```

## 本地开发（不用 Docker）

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

运行测试（需额外安装 `pytest` 与 `httpx`）：

```bash
python3 -m pytest tests/ -q
```

## 项目结构

```
app/
  main.py      # FastAPI 应用与路由（/health、/digest、/partial-digest、/design-panel）
  digest.py    # 核心酶切逻辑：IUPAC 匹配、切点计算、片段生成
  partial.py   # 部分酶切：有效切点合并、片段长度期望的精确 O(K²) 计算
  panel.py     # 面板设计：单酶切条带、可区分性判定、最小酶组选择
  schemas.py   # 请求/响应模型与校验
tests/
  test_digest.py
  test_partial_digest.py
  test_panel.py
Dockerfile
docker-compose.yml
requirements.txt
```
