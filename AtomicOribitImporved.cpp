#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <ctime>
using namespace std;

#define PI 3.141592653589 //圆周率
#define a0 5.19*pow(10,-11) //玻尔半径
#define ORBIT p_4pz
#define FILENAME "4pz_imporved.csv"
#define BOXLEN 45
//BOXLEN:2pz~18, 3pz~27, 4pz~45
//Origin作图时坐标轴长：2pz~9E-10, 3pz~1.4E-9, 4pz~2.5E-9

default_random_engine random((unsigned)time(NULL));
uniform_real_distribution<double> ux(-BOXLEN*a0,BOXLEN*a0);
uniform_real_distribution<double> uy(-BOXLEN*a0,BOXLEN*a0);
uniform_real_distribution<double> uz(-BOXLEN*a0,BOXLEN*a0);
uniform_real_distribution<double> um(0,1);

/*
本程序生成完整的原子电子云图
随机生成x,y,z,转换成r, theta, phi
按照轨道记号给出rho/rhomax函数的值
与随机M值比对，若大于等于则保留，否则舍去
输出直角坐标Coordinate文件,统计数据个数
*/

struct Coordinate//记录坐标，abc对应球极坐标的r，theta，phi或对应直角坐标的x，y，z
{
    double a = 0, b = 0, c = 0;
};

Coordinate sph2cart(Coordinate c)//球极坐标转换为直角坐标
{
    Coordinate tempc;
        tempc.a = c.a*sin(c.b)*cos(c.c);
        tempc.b = c.a*sin(c.b)*sin(c.c);
        tempc.c = c.a*cos(c.b);
    return tempc;
}

Coordinate cart2sph(Coordinate& c)
{
    double r = sqrt(pow(c.a,2) + pow(c.b,2) + pow(c.c,2));
    Coordinate tempc;
        tempc.a = r;
        tempc.b = acos(c.c/r);
        if(c.a > 0)
        {
            if(c.b > 0)
                tempc.c = atan(c.b/c.a);
            else
                tempc.c = atan(c.b/c.a) + 2*PI;
        }
        else
            tempc.c = atan(c.b/c.a) + PI;
    return tempc;
}

ostream& operator<<(ostream& out, Coordinate& c)//重载<<，直接输出结构体Coordinate
{
    out<<c.a<<','<<c.b<<','<<c.c;
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
    double rho = c.a/(2*a0);
    return pow((pow(rho,3)-10*pow(rho,2)+20*pow(rho,1)),2)*exp(-rho)*pow(cos(c.b),2)/(20480*PI*pow(a0,3))
        /(pow((pow(0.8485,3)-10*pow(0.8485,2)+20*pow(0.8485,1)),2)*exp(-0.8485)/(20480*PI*pow(a0,3)));
}

bool Criterion(Coordinate& c)//判断随机点是否符合要求
{
    double p = ORBIT(c);
    double M = um(random);
    if(p >= M)
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
    do
    {
        Coordinate sph_c, cart_c;
            cart_c.a = ux(random);
            cart_c.b = uy(random);
            cart_c.c = uz(random);
        sph_c = cart2sph(cart_c);
        if(Criterion(sph_c))
        {
            count++;
            cart_c = sph2cart(sph_c);
            outfile<<cart_c<<endl;
        }
    }while(count<100000);
    outfile.close();
    return 0;
}