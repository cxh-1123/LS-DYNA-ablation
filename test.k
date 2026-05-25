$ ============================================================
$  test.k  -  LS-DYNA 最小环境验证模型
$  说明：单个六面体单元 + 弹性材料 + 底面 4 节点固定
$        仅用于验证 lsdyna 求解器能否正常启动并跑到 normal termination。
$  单位制：mm-ms-kg-kN
$        密度 7.85e-6 kg/mm^3
$        E    = 2.10e+2 kN/mm^2 (= 210 GPa)
$ ============================================================
*KEYWORD
*TITLE
Minimal LS-DYNA Environment Test
$
$ -------- 时间控制 --------
*CONTROL_TERMINATION
$#  endtim
     1.0E-3
$
*CONTROL_TIMESTEP
$#  dtinit    tssfac
       0.0       0.9
$
*DATABASE_BINARY_D3PLOT
$#      dt
    1.0E-4
$
$ -------- Part / Section / Material --------
*PART
$#                                                                       title
Minimal Test Part
$#     pid     secid       mid
         1         1         1
$
*SECTION_SOLID
$#   secid    elform
         1         1
$
*MAT_ELASTIC
$#     mid        ro         e        pr
         1  7.850E-6  2.100E+2      0.30
$
$ -------- 几何（单位立方体）--------
*NODE
$#   nid               x               y               z
       1             0.0             0.0             0.0
       2             1.0             0.0             0.0
       3             1.0             1.0             0.0
       4             0.0             1.0             0.0
       5             0.0             0.0             1.0
       6             1.0             0.0             1.0
       7             1.0             1.0             1.0
       8             0.0             1.0             1.0
$
*ELEMENT_SOLID
$#   eid     pid      n1      n2      n3      n4      n5      n6      n7      n8
       1       1       1       2       3       4       5       6       7       8
$
$ -------- 边界条件：底面 4 节点固定 --------
*BOUNDARY_SPC_NODE
$#   nid     cid    dofx    dofy    dofz   dofrx   dofry   dofrz
       1       0       1       1       1       0       0       0
       2       0       0       1       1       0       0       0
       3       0       0       0       1       0       0       0
       4       0       1       0       1       0       0       0
$
*END
