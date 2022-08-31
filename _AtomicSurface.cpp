#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#define PI 3.141592653589 //圆周率
#define a0 5.19*pow(10,-11) //玻尔半径
#define RESOLUTION 300
#define ORBIT p_4pz
#define FILENAME "4pz_surface.csv"

/*
本程序生成轨道轮廓图
在三维空间x,y,z等距取点，转化为r,theta,phi，计算psi^2值
取20±1%等值点，转换回x,y,z，输出绘图
*/

struct Coordinate//记录坐标，abc对应球极坐标的r，theta，phi或对应直角坐标的x，y，z
{
    double a = 0, b = 0, c = 0;
};

Coordinate sph2cart(Coordinate& c)//球极坐标转换为直角坐标
{
    Coordinate tempc;
        tempc.a = c.a*sin(c.b)*cos(c.c);
        tempc.b = c.a*sin(c.b)*sin(c.c);
        tempc.c = c.a*cos(c.b);
    return tempc;
}

Coordinate cart2sph(Coordinate& c)
{
    Coordinate tempc;
        tempc.a = sqrt(pow(c.a,2) + pow(c.b,2) + pow(c.c,2));
        tempc.b = atan(sqrt(pow(c.a,2) + pow(c.b,2))/c.c);
        tempc.c = atan(c.b/c.a);
    return tempc;
}

ostream& operator<<(ostream& out, Coordinate& c)//重载<<，直接输出结构体Coordinate
{
    out<<c.a<<' '<<c.b<<' '<<c.c;
    return out;
}

double p_2pz(Coordinate& c)// psi^2(210)/psi^2(210)max
{
    double rho = c.a/(a0);//务必打括号，注意此处是宏替换
    return (pow(rho,2)*exp(-rho)*pow(cos(c.b),2)/(32*PI*pow(a0,3)))
        /(pow(2,2)*exp(-2)/(32*PI*pow(a0,3)));
}

double p_3pz(Coordinate& c)// psi^2(310)/psi^2(310)max
{
    double rho = 2*c.a/(3*a0);
    return (pow(rho,4)-8*pow(rho,3)+16*pow(rho,2))*exp(-rho)*pow(cos(c.b),2)/(648*PI*pow(a0,3))
        /((pow(1.17157,4)-8*pow(1.17157,3)+16*pow(1.17157,2))*exp(-1.17157)/(648*PI*pow(a0,3)));
}

double p_4pz(Coordinate& c)// psi^2(410)/psi^2(410)max
{
    double rho = 2*c.a/(a0);
    return pow((pow(rho,3)-10*pow(rho,2)+20*pow(rho,1)),2)*exp(-rho)*pow(cos(c.b),2)/(20480*PI*pow(a0,3))
        /(pow((pow(1.6969,3)-10*pow(1.6969,2)+20*pow(1.6969,1)),2)*exp(-1.6969)/(20480*PI*pow(a0,3)));
}

bool Criterion(Coordinate& c)//判断点是否符合要求
{
    double p = ORBIT(c);
    if(p <= 0.21 && p >= 0.14)
        return true;
    else
        return false;
}

int main()
{
    //初始化
    int count=0;
    fstream outfile(FILENAME, ios::out |ios::trunc);
    if(!outfile.is_open())
        cout<<"错误：未能成功打开文件"<<endl;

    Coordinate cart_c, sph_c;
    for(int i=0; i<RESOLUTION; i++)
    for(int j=0; j<RESOLUTION; j++)
    for(int k=0; k<RESOLUTION; k++)
    {
        cart_c.a = i/double(RESOLUTION)*18*a0-9*a0;
        cart_c.b = j/double(RESOLUTION)*18*a0-9*a0;
        cart_c.c = k/double(RESOLUTION)*18*a0-9*a0;
        sph_c = cart2sph(cart_c);
        if(Criterion(sph_c))
            outfile<<cart_c<<endl;
    }
    outfile.close();
    return 0;
}