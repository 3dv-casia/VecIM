1.依赖库：
  CGAL 4.11+
  opencv
  gmp
  glog

2.第三方库
  gco-v3.0
  scip

3.程序编译
  cd build
  cmake  -DCMAKE_BUILD_TYPE=Release ..
  cd Facade
  make
  
3.程序执行
  ./Facade    xxx/xxx/xxx/
  可执行文件  数据目录（目录下存放所有的输入输出数据）
  数据目录下需要有配置文件"config_modeling.xml"，文件内控制程序的执行流程和参数。
