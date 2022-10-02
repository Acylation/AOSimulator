# AOSimulator
 Monte Carlo simulation of atom orbits based on Schrodinger Equation

---
### Schrodinger Equation
### Monte Carlo Simulation
### Random Engine Rather than Rand()
### Todo: UI Design & Tuning
### Todo: Data Visualization by Qwt or MathGL
### Todo: Speed Up & Animation
---
### 电子云黑点图和二维切片的作图说明  

我们使用蒙特卡洛法生成电子云黑点图  

随机生成$x, y, z$三个参量，转换成极坐标，代入$ψ^2(r, θ, Φ)$中，求得$ψ^2/ψ^2_{max}$  

将$ψ^2/ψ^2_{max}$与一个0~1的随机数M作比较。若$ψ^2/ψ^2_{max}>M$，则该点有效，保留坐标。空间中每个点被保留的概率（$ψ^2/ψ^2_{max}>M$的概率）恰是$ψ^2/ψ^2_{max}$  

这样一来，就用黑点的疏密形象地表示出了$ψ^2/ψ^2_{max}$值的大小  

将生成的坐标文件导入到Origin中，作三维散点图，就能得到对应轨道的电子云黑点图  

在前述方法的基础上，筛选出靠近原点的一层平面，即得到电子云切片，如，$2p_z$轨道的$xz$切片，就是取$y$值在$-0.5a_0$ ~ $0.5a_0$之间的点，之后以$x$为横坐标，$z$为纵坐标作图得到的。

轨道轮廓图：等$|ψ|$值图，筛选上面生成数据中，$|ψ|$值在给定值附近的点形成球壳轮廓

界面图：该轮廓内电子出现的概率$p$是一特定值。界面图可以采用蒙特卡洛积分法绘出。由$n$个点形成的黑点图中，由于界面内任意点的$ψ^2$值必大于该界面的$ψ^2$值，因此只需保留$ψ^2$值最大的$n×p$个点，我们就能得到$p$对应的临界$ψ^2$值，进而参考上一段的方法作出电子云界面图

Latex Example  
$2p_z: θ=arccos(\frac{4\sqrt{2πa^5}}rψ_0e^\frac{r}{2a})$
