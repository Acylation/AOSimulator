#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
using namespace std;

#define PI 3.141592653589 //圆周率
#define a0 5.19*pow(10,-11) //玻尔半径
#define ORBIT p_3pz
#define FILENAME "3pz.csv"

/*
本程序生成完整的原子电子云图

随机生成r,theta,phi
按照轨道记号给出rho/rhomax函数的值
与随机M值比对，若大于等于则保留，否则舍去
转换极坐标与直角坐标
输出Coordinate文件,统计数据个数
*/

struct Coordinate//记录坐标，abc对应球极坐标的r，theta，phi或对应直角坐标的x，y，z
{
    double a = 0, b = 0, c = 0;
};

ostream& operator<<(ostream& out, Coordinate& c)//重载<<，直接输出结构体Coordinate
{
    out<<c.a<<' '<<c.b<<' '<<c.c;
    return out;
}

double random()//返回一个0~1之间的double数
{
    double t;
    t = rand()/double(RAND_MAX);//系统的rand()返回一个0~RAND_MAX间的整数
    return t;
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

bool Criterion(Coordinate& c)//判断随机点是否符合要求
{
    double p = ORBIT(c);
    double M = random();
    if(p >= M)
        return true;
    else
        return false;
}

Coordinate sph2cart(Coordinate c)//球极坐标转换为直角坐标
{
    Coordinate tempc;
        tempc.a = c.a*sin(c.b)*cos(c.c);
        tempc.b = c.a*sin(c.b)*sin(c.c);
        tempc.c = c.a*cos(c.b);
    return tempc;
}

int main()
{
    //初始化
    int count=0;
    srand((unsigned)time(NULL));
    fstream outfile(FILENAME, ios::out |ios::trunc);
    if(!outfile.is_open())
        cout<<"错误：未能成功打开文件"<<endl;

    do
    {
        Coordinate sph_c, cart_c;
            sph_c.a = 9*a0*random();
            sph_c.b = PI*random();
            sph_c.c = 2*PI*random();
        if(Criterion(sph_c))
        {
            count++;
            cart_c = sph2cart(sph_c);
            outfile<<cart_c<<endl;
            //outfile<<count<<' '<<cart_c<<endl;
        }
    }while(count<100000);
    outfile.close();
    return 0;
}