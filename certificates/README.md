# 全部使用免费软件的必要计算附件

本版替换此前包含 Magma 待运行项的附件。运行只需免费的 SageMath；所需 Python、PARI/GP 等由 Sage 环境提供，不需要购买 Magma，也不需要在线计算器或合作者提供授权。

采用了替代证明和较小的充分计算，而不是把 Magma 命令机械翻译成 Sage。所有保留的计算均已实际执行。正文仍提供全局界、系数引理和引用定理；本附件不是全文的形式化证明。

## 原来四组 Magma 任务如何解决

| 原任务 | 本版处理 | 需要运行的文件 |
|---|---|---|
| D=21 的空间维数和三个 Hecke 行列式 | 比较指数 3、4、5 的系数，直接排除所有源权重；不再需要维数或行列式 | 无；完整证明见 notes/d21_local_free.tex |
| D=169、361 的权重二端点 | 精确 CM 类数、Shimizu 面积和亏格公式给出 S2=0 | cubic_weight_two_genus.sage |
| D=725 的权重二端点 | Eichler 类数公式给出四元数最大序类数 1，从而 S2=0 | quartic_725_weight_two.sage |
| D=1125 的剩余端点 | 证明 S2=0，并构造至少 5 维的 S5(epsilon) 子空间；足以排除剩下的权重三来源 | quartic_1125_sufficient.sage |
| D=12 的两个维数 | 已发表环结构、明确的低权单项式和常数项计算给出 S4(1)=2、S5(epsilon)=1 | 无；完整证明见 notes/d12_ring_dimensions.md |

尤其注意：本版没有声称重新算出了原稿的三个行列式，也没有声称验证了 D=1125 的精确维数 S3=2、S5=6。新证明不再需要这些数值。D=1125 应将原来的精确维数陈述改为本版已经证明的充分结论。

## 保留的原算术计算

- e3_cusp_check.sage：E3 尖点分支的数域、本原特征和特殊值。
- e2_square_branch_check.sage：二进平方分支的完整枚举及不等式。
- e2_high_degree_check.sage：七、八次域的有限检查；七次域独立枚举，八次域名单完备性仍使用正文引文。
- e2_candidate_enumeration.sage：完整枚举 771 个域，保留所有筛选记录。
- e2_allclass_balanced.sage 与 e2_indecomposable_search.sage：任意窄类分量的格证书，共 15 个域；D=3969 另存完整 18 点检查。
- e2_dyadic_bound.sage：特殊值、正分解和权重截断。
- rq_certificate.py 与 rq_e2_exact.sage：实二次域完整缩减和局部关系。

根目录 JSON 和 logs/ 是实际输出；notes/ 是数学依据。格搜索 JSON 中的 unresolved_search 仅表示该有界格搜索没有取得证书，不表示合并后的证明仍有缺口；其三个域已由新增免费端点证书处理。

## 运行方法

使用 SageMath 10.8 或兼容环境。本次实际计算使用 passagemath-standard 10.8.12。
SageMath 官方安装说明：
https://doc.sagemath.org/html/en/installation/index.html

解压后，在终端进入 free_certificates 文件夹。只复算此次新增的三个端点计算：

    sage -python run_all.py --endpoints-only

完整复算所有保留程序：

    sage -python run_all.py

新输出写到 rerun/，不会覆盖附带的原始输出。运行成功后分别显示 FREE_ENDPOINTS_PASSED 或 FREE_ARITHMETIC_SUITE_PASSED，并生成对应状态 JSON。

例如在 Windows 中，把文件夹放到 C:\math\free_certificates；打开已安装 Sage 的 Ubuntu/WSL 后执行：

    conda activate sage
    cd /mnt/c/math/free_certificates
    sage -python run_all.py --endpoints-only

这些命令在 Ubuntu 终端执行，不是在 Sage 的交互提示符或网页输入框内执行。

## 正文要同步修改

采用 manuscript_changes.md 中的修改清单，尤其替换旧 Magma 计算断言。D=21 和 D=12 的计算已经完整写在证明中，所以按“正文算清就不保留重复程序”的要求，不再附上用于内部核算的两个短程序。

新证明引用 Shimizu/Voight 的亏格公式、Kirschmer–Voight 的类数公式、Jacquet–Langlands 对应，以及 Aoki/Dursthoff 的环结构。完整出处与适用条件列在相应 notes/ 文件中；免费复算这些有限算术并不替代对引用定理的使用。
