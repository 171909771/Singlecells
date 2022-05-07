# 安装velocyto.R
```r
install_github("velocyto-team/velocyto.R")
```
## 当出现 “/usr/bin/ld: cannot find -lboost_filesystem” 和 “/usr/bin/ld: cannot find -lboost_system” 问题时
### 有的问题 apt-get 解决不了，必须使用 aptitude 解决，有的问题，用 aptitude 解决不了，必须使用 apt-get
```linux
sudo aptitude install libboost-all-dev
```
