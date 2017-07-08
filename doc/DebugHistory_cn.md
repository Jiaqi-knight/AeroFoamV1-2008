# AeroFoamV1-2008移植

这个似乎会简单一些。

## 建立repository

```shell
wget https://home.aero.polimi.it/freecase/?download=AeroFoam31072008.tar.gz -O AeroFoam31072008.tar.gz 
mkdir AeroFoamV1-2008
cd AeroFoamV1-2008
mkdir src
cd src
tar xf ../../AeroFoam31072008.tar.gz
mv AeroFoam/* .
rmdir AeroFoam
cd ..
mkdir tutorials
mv src/Data\ of\ CAE3D\ Example\ AGARD445.6/ tutorials/
wget https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-dev/master/.gitignore -O .gitignore
git add .gitignore
git add *
git status
git commit -m "init"
git remote add origin https://github.com/chengdi123000/AeroFoamV1-2008.git
git remote -v
git config credential.helper store #记住用户名密码，仅针对repository这个有效。
git push -u origin master
```

## 编译调试

- 添加readme.md文件，指明这个项目的来源和目的，并声明自己不是原作者。

- 修改Make/files和Make/options，使其能够编译。

- ```shell
  grep -n "\(boundaryField()\)\(.*=.\+\)" *.* #查看所有boundaryField()后面有等于符号的地方，说明对边界场进行了赋值，但是边界场赋值在openfoam 4.1中必须以boundaryFieldRef()的形式出现。原来的API不区分引用和赋值的做法被废弃了。
  # 然后可以用命令sed进行批量替换
  sed -i -e "s/\(boundaryField()\)\(.*=.\+\)/boundaryFieldRef()\\2/g" *.*
  ## 上述表达式`s/原有字符串模式/新字符串模式/g`的意思是
  ## s: 替换
  ## g: greed，贪婪模式，尽可能匹配得多
  ## 原有字符串模式，分为两段，前面是`\(boundaryField()\)`，后面是`\(.*-.\+\)`
  ### 这里有大量转义字符`\`，去掉转义字符是正常的字符串："s/(boundaryField())(.*=.+)/boundaryFieldRef()\2/g"
  ### 第一段的意思是匹配boundaryField()yyy=xxx，用括号括起来表示这是一个独立的pattern，与后面的pattern相组合。在最后的新字符串模式中，这个pattern可以用\1来引用
  ### 第二段的意思是：`.`匹配任意字符，`.*`表示任意字符串重复任意多次，然后是有个等于符号`=`，后面`.\+`的意思
  # 反正比较恶心的是这里有两层转义，一层是shell的转义，一层是sed所用regular expression的转义。
  ```

- wmake之后发现还是有很多error，但只集中在两个文件里：`CreateConnectivity.C`中和`CAE3DTools.C`

  - `CreateConnectivity.C`的问题是fvMesh中没有allFaces(), allNeighbour(), allOwner()函数。

  - ```shell
    di@ubuntu:~/OpenFOAM/di-4.1/solvers/AeroFoamV1-2008/src$ grep allFaces CreateConnectivity.C
    labelField id_LL_i( mesh.allFaces().size(), 0 );
    labelField id_RR_i( mesh.allFaces().size(), 0 );
    id_adv = ( mesh.allFaces().size() - mesh.Sf().size() - 1 )/10;
                printf("%3.0f %% \n", ( i - mesh.Sf().size() - 1 )*100.0/( mesh.allFaces().size() - mesh.Sf().size() - 1 ) );
                id_adv = id_adv + ( mesh.allFaces().size() - mesh.Sf().size() )/10;
    di@ubuntu:~/OpenFOAM/di-4.1/solvers/AeroFoamV1-2008/src$ grep allOwner CreateConnectivity.C
            id_L  = mesh.allOwner()[i];
            id_L  = mesh.allOwner()[i];
    di@ubuntu:~/OpenFOAM/di-4.1/solvers/AeroFoamV1-2008/src$ grep allNeighbour CreateConnectivity.C
            id_R  = mesh.allNeighbour()[i];
            id_R  = mesh.allNeighbour()[i];
    ```

  - 可见allFaces()是用于取得size的，allNeighbour()和allOwner()是用于取得一些id_L, id_R的。

  - 解决方案

    - 查找OpenFOAM 1.5-dev的代码，没找到，用[openfoam-extend-Core-OpenFOAM-1.5-dev](https://github.com/Unofficial-Extend-Project-Mirror/openfoam-extend-Core-OpenFOAM-1.5-dev)替代一下，发现allFaces()的完全定义是：

    - ```c++
      class polyMesh
      :
          public objectRegistry,
          public primitiveMesh
      {
                  //- Return all faces, including inactive ones
                  const faceList& allFaces() const;
      }
      ```

    - 从这个OpenFOAM 1.1版本的functions来看（[list of all class members](http://openfoamcfd.sourceforge.net/doc/Doxygen/html/functions.html) ）来看，在老版本中：

      - faceOwner()faceNeighbour()是primitiveMesh类实现的
      - allOwner()/allNeighbour()是polyMesh类实现的
      - owner()/neighbour()是fvMesh类实现的。
      - allFaces()是polyMesh实现的，polyMesh还额外实现了一个faces()函数，似乎faces()是allFaces()中的一部分。

    - 通过比较of41和of15-extend的代码，特别是polyMesh的构造函数，可以发现，of15中的allFaces()就是of41中的faces()，所以mesh.allFaces().size()的语义是所有的（从`constant/polyMesh/faces`文件中读取到的）面的数量，所以可以替换mesh.allFaces().size() 为mesh.faces().size()；

    - 但是of15-extend中已经没有了allOwner()，说明这份AeroFoam的代码比我们想象的还要老！还得去找更老的代码才能知道这是什么意思。

    - 但是代码很难找，找到了（https://sourceforge.net/projects/foam/files/foam/）还得解压出来看，我想先改一下再说。把allOwner()直接改成faceOwner()试试，因为这是直接返回内部owner_变量的函数。

    - 最终的sed代码如下：

    - ```shell
      sed -i -e "s/mesh.allFaces().size()/mesh.faces().size()/g"  CreateConnectivity.C
      sed -i -e "s/mesh.allOwner()/mesh.faceOwner()/g"  CreateConnectivity.C
      sed -i -e "s/mesh.allNeighbour()/mesh.faceNeighbour()/g"  CreateConnectivity.C
      wmake #至少形式上没有CreateConnectivity.C的错误了。
      ```

  - 然后再来看`CAE3DTools.C`的问题，打开这个文件就看见这么一句话，一口老血喷出三丈远：

  - ```c++
        // TODO: *** Read structural mesh data from .mail file of Code_Aster *** [???]

        // Read StructuralModel.vertices and initialize memory
        f1 = fopen("./Data/StructuralModel.vertices", "r");
        if ( f1 == NULL ) Info << "ERROR: File StructuralModel.vertices not found!" << nl;
        fscanf(f1, "%i", &N);
        (*xx_s_v) = vectorField(N);

    ```

  - 说明两个问题：

    - 这个版本与Code Aster之间的连接还不正常。
    - 结构模型数据的路径是写死的。

  - 然后定位到错误，这时我想到一个问题，为什么`CAE3DTools.C`有问题而`CAE2DTools.C`没有问题。然后我发现2D版本没有用Matrix类。看来都是Matrix类惹的鬼。

  - 然后统计一下，42个错误中有15个和`Matrix<scalar>`类型有关。

    - 查阅of41的文档，of41中的`Matrix<>`的接口应该是`Matrix<Form<Type>,Type>`，改！

    - ```shell
      grep -n "Matrix<scalar>"  CAE3DTools.C
      sed -i -e "s/Matrix<scalar>/Matrix<SquareMatrix<scalar>,scalar>/g" CAE3DTools.C
      wmake|less
      ```

    - OK！

  - 然后看看`Utilities`，这是后处理工具，主要是做图和显示位移用的。

    - showDisplacement是一个单独的小程序，改改Make/files的exe生成的位置，改改boundaryField就好了。
    - gnuplot有两个版本的，都挺老的了，用`apt-cache show gnuplot-qt`，ubuntu 16.04自带的gnuplot已经到了4.6.6-3。
    - 远程ssh的时候需要有才能有X11映射才能显示图形，在windows下建议选择`vcxsrv`，链接https://sourceforge.net/projects/vcxsrv/，比xming要好一些。
    - 运行了一下相关的脚本，发现至少需要一些算例才行，那就先放一边吧，先下载一个算例来看看。

## 运行

首先去下载一些算例，在freecase的网站上有不少呢。

先来个最简单的斜激波吧。

```shell
# from https://home.aero.polimi.it/freecase/?OpenFOAM_%2B_Code_Aster:Aerodynamic_problems:Oblique_shock
wget "https://home.aero.polimi.it/freecase/?download=ObliqueShock.tar.gz" -O ObliqueShock.tar.gz
tar xf ObliqueShock.tar.gz
# 解压出来有 ObliqueShock_AeroFoam  ObliqueShock_rhopSonicFoam  ObliqueShock_rhoSonicFoam  ObliqueShock_sonicFoam  四个文件夹
cd ObliqueShock_AeroFoam
# 有从细到粗多个网格的算例
# ObliqueShock_A  ObliqueShock_B  ObliqueShock_C  ObliqueShock_D
cd ObliqueShock_A
tree # 查看文件树
#.
#├── Data
#│   └── Readme #Data文件夹存放connectivity文件，如果有气动弹性计算，也存放气弹模型和数据
#├── Log
#│   └── Readme #Log文件夹存放Log，可以用utilities中的工具做图。
#└── ObliqueShock_A
#    ├── 0
#    │   ├── p
#    │   ├── T
#    │   └── U
#    ├── 0.001
#    │   ├── eet_tilde.gz
#    │   ├── mm.gz
#    │   ├── p.gz
#    │   ├── rrho.gz
#    │   ├── T.gz
#    │   ├── U.gz
#    │   └── uniform
#    │       └── time.gz
#    ...
#    ├── constant
#    │   ├── polyMesh
#    │   │   ├── blockMeshDict #早期的blockMeshDict似乎是放这里的。
#    │   │   ├── boundary
#    │   │   ├── faces
#    │   │   ├── neighbour
#    │   │   ├── owner
#    │   │   └── points
#    │   └── thermodynamicProperties
#    └── system
#        ├── controlDict
#        ├── fvSchemes #可以看下这些文件。
#        └── fvSolution
```

然后运行出错了！！！信息中没有出错位置的信息，怎么办呢，加调试选项！

```shell
# src/Make/options
# 加入-g选项，开启调试
EXE_INC = -g \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude

EXE_LIBS = -lfiniteVolume \
           -lmeshTools
```

在wmake之后，运行，提示错误信息：

```shell
#0  Foam::error::printStack(Foam::Ostream&) at ??:?
#1  Foam::sigSegv::sigHandler(int) at ??:?
#2  ? in "/lib/x86_64-linux-gnu/libc.so.6"
#3  __fprintf_chk in "/lib/x86_64-linux-gnu/libc.so.6"
#4  ? at ~/OpenFOAM/di-4.1/solvers/AeroFoamV1-2008/src/AeroFoam.C:635
#5  __libc_start_main in "/lib/x86_64-linux-gnu/libc.so.6"
#6  ? at ??:?
```

然后看看AeroFoam.C的第635行是什么内容：

``` c++
      	fclose(ff);
```

居然是文件关错了。等等，看了一下目录，你也没有成功打开新的文件呀，看来真是有问题了。尝试进gdb模式调试一下。

```shell
gdb AeroFoam
(gdb) b 635
(gdb) run
# 出现错误之后，先看看栈
(gdb) bt
#0  ___fprintf_chk (fp=fp@entry=0x0, flag=flag@entry=1,
    format=format@entry=0x4803e7 "%e %e %e\n") at fprintf_chk.c:30
#1  0x000000000041d57c in fprintf (__fmt=0x4803e7 "%e %e %e\n", __stream=0x0)
    at /usr/include/x86_64-linux-gnu/bits/stdio2.h:98
#2  main (argc=1, argv=0x7fffffffd8f8) at AeroFoam.C:634
## 那么我们先移动到main函数中看看
(gdb) up
#1  0x000000000041d57c in fprintf (__fmt=0x4803e7 "%e %e %e\n", __stream=0x0)
    at /usr/include/x86_64-linux-gnu/bits/stdio2.h:98
98                              __va_arg_pack ());
(gdb) up
#2  main (argc=1, argv=0x7fffffffd8f8) at AeroFoam.C:634
634             fprintf(ff, "%e %e %e\n", err_rrho, err_mm, err_eet_tilde);
## 这里是第634行，猜测是ff==0，也就是没有能顺利打开文件去写入。
(gdb) disp ff
1: ff = (FILE *) 0x0 ##看来就是这个问题了。
## 再来看看代码
(gdb) l -2
632             //IF{ writeResidual 1
633             ff = fopen("./Log/residuals.txt", "a");
634             fprintf(ff, "%e %e %e\n", err_rrho, err_mm, err_eet_tilde);
635             fclose(ff);
636             //FI} writeResidual 1
637
638             // Update solution at time level (k-1)->oo and (k)->o
639             rrho_oo      = rrho_o;
640             mm_oo        = mm_o;
641             eet_tilde_oo = eet_tilde_o;
## fopen的"a"模式是什么鬼？
## 查找资料，"a"表示在文件末尾添加内容，如果没有文件则创建文件。这个模式本身没啥问题。
## 可能问题出在"Log"目录上，会不会是和没有Log目录有关系。
(gdb) quit ## 出去尝试一下。
```

```shell
## 将上一个目录的Log和Data文件夹挪进来！
mv ../Data ./
mv ../Log ./
foamJob AeroFoam ## foamJob是OpenFOAM标准的运行脚本
## 查看输出
tail -f log
#<CTRL>+c 退出
## 调用一下后处理工具
../../../../utilities/GNUplot4.2/plotStats
#
#Cannot open load file 'convergence.plt'
#"convergence.plt", line 0: util.c: No such file or directory
#
## 原来是因为没有找到convergence.plt呀，看来有一点儿路径设置的问题。
```

参考了一下网上的这个回答：[Getting the source directory of a Bash script from within](https://stackoverflow.com/a/59916/4592964)，把`utilities/GNUplot4.2/plotStats`给修改了一下：

```shell
#!/bin/bash

if [ $# -eq 0 ]; then # `$#` denotes number of input
todo="residual"
else
todo=$1
fi

# modification here: get the dirname of the script.
# Not where it is invoked
DIR=`dirname $0`
if [ $todo = "residual" ]; then
gnuplot $DIR/convergence.plt
fi
if [ $todo = "2D" ]; then
gnuplot $DIR/loads2D.plt 2&> /dev/null	
gnuplot $DIR/movement2D.plt 2&> /dev/null
python $DIR/xCp2D.py
fi
if [ $todo = "3D" ]; then
gnuplot $DIR/loads3D.plt 2&> /dev/null
gnuplot $DIR/modal3D.plt 2&> /dev/null
echo "Press return to export current image into deform3D.eps file..."
python $DIR/xyCp3D.py
fi
```

然后发现`convergence.plt`中居然又把`residuals.txt`的路径写死了。为了尽可能的少动代码，我得把后处理工具都给改一下才行。

```shell
cd utilities
rm -rf GNUplot4.0
mv GNUplot4.2 GNUplot
cd GNUplot
grep "\.\.\/Log" *.*
sed -i -e  "s/\.\.\/Log/\.\/Log/g" *.*
cd ../showDisplacement
```

并且把文件夹也移动了一下，成这样：

```shell
cd ObliqueShock_A
mv ObliqueShock_A/* .
rmdir ObliqueShock_A/
mkdir Fig
cd ..

cd ObliqueShock_B
mv ObliqueShock_B/* .
rmdir ObliqueShock_B/
mkdir Fig
cd ..

cd ObliqueShock_C
mv ObliqueShock_C/* .
rmdir ObliqueShock_C/
mkdir Fig
cd ..

cd ObliqueShock_D
mv ObliqueShock_D/* .
rmdir ObliqueShock_D/
mkdir Fig
cd ..
```

我觉得一个很重要的原因是，c的api函数fopen只能在已经存在的文件夹下建立文件，不能直接新建路径中的文件夹，如果fopen("./dir/file")而"./dir"还不存在，就会报错！

用`foamJob AeroFoam`运行，但是只有ObliqueShock_D运行成功。其他均出现和rrho有关的段错误（segmentation fault）。这TMD真是日了狗了。应该是数组越界了。

但是gdb调试发现数组的index被优化掉了，还得重新编译，加入`-O0`选项，覆盖掉`-O3`选项（后置的`-O`选项可以覆盖前面的）。再次进入gdb，可以看出变量确实越界了：

```shell
(gdb) l
230                     n_i_        = n_i_/mag(n_i_);
231
232                     // Find id_LL = (j - 1) ( see Jameson )
233                     id_LL        = (*id_LL_i)[i];
234                     id_RR        = id_R;
235                     rho_LL_      = (*rrho)[id_LL];
236                     m_LL_        = (*mm)[id_LL];
237                     et_tilde_LL_ = (*eet_tilde)[id_LL];
238
239                     // Cells CGs
(gdb) print i
$1 = -750
(gdb) print ii
$2 = 0
```

根据代码中的说明，`ii`表示局部index，`i`表示全局index，但是明显`whichFace`的语义是输入全局index，返回局部index，这里应该是弄反了。

进一步地，通过搜索代码可以发现，OpenFOAM 2.x的时候求解器中还会使用`whichFace()`函数，到OF41的时候就没有自带的求解器使用这个函数了，说明这个函数其实是比较老的函数了。

再进一步看`forAll`宏的定义，这玩意是从0开始的。所以whichFace(0)会返回负值。

```c++
define forAll(list, i) \
  for (Foam::label i=0, i<(list).size(); i++)
```

于是尝试用以下方式来解决：

```shell
#查看一下所有用了whichFace的地方
grep "whichFace(ii)" *.* 
#试验一下sed语句 
echo "ia   = mesh.boundaryMesh()[iPatch].whichFace(ii);"|sed -e "s/whichFace(ii)/start()+ii/g"
# 替换所有whichFace(ii)为start()+ii
sed -i -e "s/whichFace(ii)/start()+ii/g" *.*
```

后来我才发现[AeroFoam QuickStart](https://home.aero.polimi.it/freecase/?download=AeroFoam_QuickStart.pdf)中首先就要求改`whichFace`这个函数的定义。但是我感觉作者是弄错了，这是他们程序的bug，而不是OpenFOAM程序包的bug。

修改之后`wmake`即可，但是注意把加入的`-O0`选项去掉，不然运行巨慢！

收敛残差图：

![收敛残差图](./imgs/AeroFoam_ObliqueShock_A_residuals.png)

## 算例

### ObliqueShock

同上，下载运行即可

###  IncompressibleLimit

```shell
wget "https://home.aero.polimi.it/freecase/?download=IncompressibleLimit.tar.gz" -O IncompressibleLimit.tar.gz
tar xf IncompressibleLimit.tar.gz
rm IncompressibleLimit.tar.gz
cd IncompressibleLimit
mv Data Cylinder_15k_M010/
mv Log Cylinder_15k_M010/
mkdir -p Cylinder_15k_M010/Fig
git add ../IncompressibleLimit
git commit -m "IncompressibleLimit initial add"
cd Cylinder_15k_M010
foamJob AeroFoam #everything should be OK. However, the time is too long. 
```

###  NACA 0012

```shell
wget https://home.aero.polimi.it/freecase/?download=NACA0012.tar.gz -O NACA0012.tar.gz
mkdir NACA0012
mv NACA0012.tar.gz NACA0012
cd NACA0012
tar xf NACA0012.tar.gz
rm NACA0012.tar.gz
```

运行似乎也没有什么问题

### 二维薄翼型

```shell
wget https://home.aero.polimi.it/freecase/?download=ThinAirfoil.tar.gz -O ThinAirfoil.tar.gz
```

同理，运行正常！

可以采用如下方式运行：

```shell
foamJob -case ThinAirfoil_7k_M120_a1 AeroFoam
# 这样可以避免不断地cd操作。
```

### OneraM6

```shell
mkdir OneraM6
cd OneraM6
wget https://home.aero.polimi.it/freecase/?download=OneraM6.tar.gz -O OneraM6.tar.gz
tar xf OneraM6.tar.gz
rm OneraM6.tar.gz
mv OneraM6/* .
rmdir OneraM6
mv OneraM6_340k_M084_a3/* .
rmdir OneraM6_340k_M084_a3/
git add .
foamJob AeroFoam 
```

应该注意到，这个算例设置的`loadConnectivity=1`，所以Data文件夹下必须存在`LR2LLRR.txt`文件。

### Woodward & Colella算例，即前台阶

```shell
wget https://home.aero.polimi.it/freecase/?download=WoodwardColella.tar.gz -O WoodwardColella.tar.gz
tar xf WoodwardColella.tar.gz
rm WoodwardColella.tar.gz
cd Woodward\&Colella/
mv Woodward\&Colella/* .
rmdir Woodward\&Colella/
rm Data/Readme Log/Readme
git add .
git commit -m "forwardStep case added"
foamJob AeroFoam
```

### RAEplane

```shell
wget https://home.aero.polimi.it/freecase/?download=RAEplane.tar.gz -O RAEplane.tar.gz
tar xf RAEplane.tar.gz
cd RAEplane/
mv RAEplane_460k_M090_a1/* .
rmdir RAEplane_460k_M090_a1/
rm Data/Readme Log/Readme
git add .
git commit -m "RAEplane case added"
foamJob AeroFoam
```

### YF-17

```shell
wget https://home.aero.polimi.it/freecase/?download=YF17.tar.gz -O YF17.tar.gz
tar xf YF17.tar.gz
cd YF17_530k_M028_a19
mv YF17_530k_M028_a19/* .
rm Data/Readme Log/Readme
git add .
git commit -m "YF-17 case added"
foamJob AeroFoam
```

### ZPG平板边界层

```shell
wget https://home.aero.polimi.it/freecase/?download=Blasius.tar.gz -O Blasius.tar.gz
tar xf Blasius.tar.gz
cd Blasius
mv Blasius_3k_M030/* .
rmdir Blasius_3k_M030/
rm Log/Readme
git add .
git commit -m "ZPG case added"
foamJob AeroFoam
```

###  PitchingAirfoil 2D气弹问题

```shell
 wget https://home.aero.polimi.it/freecase/?download=PitchingAirfoil.tar.gz -O PitchingAirfoil.tar.gz
```

### AileronBuzz 2D气弹算例

```shell
wget https://home.aero.polimi.it/freecase/?download=AileronBuzz.tar.gz -O AileronBuzz.tar.gz
```

### AGARDWing 3D气弹算例

```shell
wget https://home.aero.polimi.it/freecase/?download=AGARDWing.tar.gz -O AGARDWing.tar.gz
tar xf AGARDWing.tar.gz
cd AGARDWing
rm ../AGARDWing.tar.gz
mv AGARDWing_M0960_q1/* .
rm AGARDWing_M0960_q1/ -r
rm Data/Readme Log/Readme
git add .
git commit -m "AGARDWing 3D aeroelasticity case added"
foamJob AeroFoam
```



# 参考文献

主要资源来自意大利米兰理工的航空科学与技术系的一些工作。

联系人：

- Giulio Romanelli: giulio.romanelli@gmail.com
- Elisa Serioli: elisa.serioli@gmail.com

## AeroFoam

1. [FreeCASE - Free(dom) Computational AeroServoElasticity](https://home.aero.polimi.it/freecase/?FreeCASE)
2. [AeroFoam QuickStart](https://home.aero.polimi.it/freecase/?download=AeroFoam_QuickStart.pdf)
3. G. Romanelli and E. Serioli. A "free" approach to the modern Computational Aeroelasticity. Politecnico di Milano. July 23, 2008. [the Presentation](https://home.aero.polimi.it/freecase/?download=RomanelliSerioli_MSCThesis_Presentation.pdf)
4. G. Romanelli, E. Serioli and P. Mantegazza. A new accurate compressible inviscid solver for aerodynamic applications. 3rd OpenFOAM Workshop by Dipartimento di Energetica, Politecnico di Milano, Milan, Italy. July 9-11, 2008. [the Presentation](https://home.aero.polimi.it/freecase/?download=RomanelliSerioliMantegazza_3OFWorkshop_Presentation.pdf)
5. 程序：AeroFoam-v2.tar.gz complete source code by G. Romanelli. Updated to January 1, 2013. [Download.](https://www.aero.polimi.it/freecase/?download=AeroFoam01012013.tar.gz) 注：sourceforge上的AeroFoam文件比这个版本老。
6. 算例：https://home.aero.polimi.it/freecase/?OpenFOAM_%2B_Code_Aster:Aerodynamic_problems

## Code Aster

1. 官网：http://www.code-aster.org/
2. 简介presentation：http://www.code-aster.org/UPLOAD/DOC/Presentation/plaquette_aster_en.pdf
3. 算例：https://home.aero.polimi.it/freecase/?OpenFOAM_%2B_Code_Aster:Structural_problems

##  Aeroservoelastic problems

1. 算例：https://home.aero.polimi.it/freecase/?OpenFOAM_%2B_Code_Aster:Aeroservoelastic_problems

## 其他

1. 实时程序接口 (Real Time Application Interface, RTAI): https://www.rtai.org/
2. 多体动力学分析软件MBDyn (MultiBody Dynamics analysis software, MBDyn): https://www.mbdyn.org/
