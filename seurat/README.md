其中涉及`vars.to.regress`的用法
-----
- https://satijalab.org/seurat/articles/cell_cycle_vignette.html # 排除细胞周期因素`vignette`

    主要涉及 `ScaleData` 及 `SCTransform` 两个函数

Cerebro shiny 可视化
---
计算cluster之间的关系，细胞周期，线粒体/核糖体占比，富集分析及GSVA，Trajectory 
- https://romanhaa.github.io/cerebroApp/articles/cerebroApp_workflow_Seurat.html

DefaultAssay
----
设置assay的默认栏

Findmarker 中的logFC取值
----
- https://zhuanlan.zhihu.com/p/569706461

公式：log2(mean(expm1(group1.x))+1)-log2(mean(expm1(group2.x))+1)

group1.x和group2.x是RNA的data数据

简而言之：log2（平均化（标准化的值））+1）
