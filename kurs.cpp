#include <cmath>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <new>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <vector>
#include <algorithm>
#include <valarray>

#include "fft.h"
using namespace std;
std::random_device __randomDevice;
std::mt19937 __randomGen(__randomDevice());
std::normal_distribution<double> __normalDistribution(1, 0.5);
typedef complex<double> w_type;
const double TwoPi = 6.283185307179586;

inline double randn()
{

    double randNum = -1;
    while (randNum < 0.0 || randNum > 1.0)
    {
        randNum = __normalDistribution(__randomGen);
    }

    return randNum;
}

void kurs()
{
    double rand[256] = { 0.3161, 0.3796, -0.3973, 0.6517, 1.0577, 0.4964, -2.3114, -0.0204, -0.1563, -0.9320, -0.6980, 0.2410, -0.3572, 0.3742, 1.0120, 0.1753, -0.2605, 1.0462, 1.7406, -0.0487, -0.2152, 0.9676, 0.8713, 0.1876, 1.1884,
    0.6226, -0.1834, 0.5574, -1.2412, 0.8256, -0.5823, 1.0304, 1.1742, -0.0580, -1.8341, -0.2718, 1.4027, -1.4113, -0.6346, -0.3051, 0.3708, 0.0384, -0.8277, 0.7659, -0.5556, 0.7656, 1.4730, 0.6280, -0.8831, 0.3749, -0.7090,
    2.2510, -0.7774, 0.4938, -0.5067, 0.0816, -0.5777, -0.1620, 0.7778, 0.6305, -0.1853, 1.0270, 0.8523, -0.9765, -0.3319, -0.8007, -0.3480, -0.1973, -0.5845, -0.2402, 0.0397, -0.5809, 1.3704, 1.8565, 0.5005, -0.5052,
    -0.6275, -0.4058, 0.7544, -0.6054, 0.3543, -0.3921, 1.5454, 0.9244, -0.4909, 1.3404, -0.8311, 1.0693, -0.2094, -0.6028, -0.9629, 0.5400, 0.5213, -1.7641,  -2.4523, 0.1134, 0.2684, -0.0339, 0.5066, -0.6224, -0.1511,
    0.6350, 0.3369, -2.1538, -0.7569, -0.6598, -0.9968,  -1.1628, 0.4475, 0.8980, -1.2692, -0.4334, 1.3587, -0.5852, 0.5043, 0.5884, 0.4292, 0.6517, 2.2047, -1.4421, -1.0515, 1.5248, -0.3696, -1.0000, 0.7456, 0.7180, -1.6072,
    -0.5354, -0.3020, -0.7623, 0.1973, -0.1311, 0.7433, -0.3598, 2.3831, -0.4824, -1.3775, -1.7841, -0.1764, 0.8817, 0.8161, 0.5286, 0.2667, 0.1775, 0.6843, -1.1475, -0.3348, 1.5890, 0.0441, -0.1317, -0.2336, 0.2670, -1.1256,
    -0.2628,-0.1233, -2.9872, 0.3039, -1.5244, -0.3453, -0.1189, 1.7669, -0.2798, 1.1211, 0.6599, 0.9366, 0.8063, -1.0746, -0.6174, 0.1416, 0.7667, 0.0856, -0.4501, 1.1193, -1.2839, 0.2175, -1.1182, -1.5717, 2.0260,
    1.6929, 1.1863, 0.9172, -0.4083, 0.9755, 1.5850, -1.6666, -1.2655, -0.6887, -1.0670, 0.0766, -0.7500, 1.7263, -1.1390, 0.2293, 0.4413, 0.5597, -0.1665, 1.3775, -0.2048, 0.0735, -0.3266, -0.6826, 0.6226, -0.1834, 0.5574, -1.2412, 0.8256, -0.5823, 1.0304, 1.1742, -0.0580, -1.8341, -0.2718, 1.4027, -1.4113, -0.6346, -0.3051, 0.3708, 0.0384, -0.8277, 0.7659, -0.5556, 0.7656, 1.4730, 0.6280, -0.8831, 0.3749, -0.7090,
    2.2510, -0.7774, 0.4938, -0.5067, 0.0816, -0.5777, -0.1620, 0.7778, 0.6305, -0.1853, 1.0270, 0.8523, -0.9765, -0.3319, -0.8007, -0.3480, -0.1973, -0.5845, -0.2402, 0.0397, -0.5809, 1.1193 };

    int maxj2 = 1000;
    int kt = 200;
    int kt_dop = kt - 2;
    int kp = 1000;
    double fi = 1.55;
    double Q = 0.03;
    double* noise;
    double* abs_err_grad = new double[maxj2];



    /*Определение fi2 методом простого перебора
    В программе генерируется  модельный зашумленный гармонический  сигнал
    с фазой fi
    производится дискретизация при количестве точек меньшей, чем количество
    периодов.Вычисляются фаза полученного вторичного сигнала
    количество периодов сигнала 100 или 102. 100 - критическое значение
    при kp = 100 вторичный сигнал близок к прямой линии, т.к.отношение
    kp / kt = 100 / 20 = 5 - число целое.Заранее выбрать "правильное" количество
    отсчетов нельзя, т.к.кол - во периодов сигнала kp неизвестно */
    for (int j2 = 0; j2 < maxj2; j2++) { //Цикл для стат.испытаний 
        //Определение количества отсчетов, не кратного количеству периодов
        noise = new double[256];
        for (int i = 0; i < kt + 1; i++)
        {
            noise[i] = randn();
            //std::cout << noise[i];

        }
        vector<w_type> y3;
        for (int i = 0; i < kt + 1; i++)
        {
            y3.push_back(sin(2 * M_PI * kp * (i - 1) / kt + fi)); //генерация модельного сигнала при кол-ве отсчетов kt (kt+1 вместо kt)  
            y3[i] += Q * noise[i]; //первый вторичный сигнал  ( kt вместо kt-1 ) 
            //std::cout << y3[i];
        }
        vector<w_type> y4;
        for (int i = 0; i < kt_dop + 1; i++) //генерация модельного сигнала при кол-ве отсчетов kt-2 (kt_dop+1 вместо kt_dop) 
        {
            y4.push_back(sin(2 * M_PI * kp * (i - 1) / kt_dop + fi)); //Второй вторичный сигнал (kt_dop вместо kt_dop+1) 
            y4[i] += Q * noise[i];
            //cout << y4[i];
        }


        auto miny3 = min_element(y3.begin(), y3.end(), [](auto a, auto b) {return a.real() < b.real(); });
        auto maxy3 = max_element(y3.begin(), y3.end(), [](auto a, auto b) {return a.real() < b.real(); });
        //cout << endl << miny3->real() << ' ' << maxy3->real() << endl;
        auto miny4 = min_element(y4.begin(), y4.end(), [](auto a, auto b) {return a.real() < b.real(); });
        auto maxy4 = max_element(y4.begin(), y4.end(), [](auto a, auto b) {return a.real() < b.real(); });
        //cout << endl << miny4->real() << ' ' << maxy4->real() << endl;
        if ((maxy3 - miny3) < (maxy4 - miny4))
        {
            kt = kt_dop; /*если диапазон амплитуд вторичного сигнала при количестве
                        отсчетов kt меньше, чем при количестве отсчетов kt - 2, то дальше будет
                        обрабатываться результат регистрации вторым АЦП, работающим с немного
                        более высокой частотой, так что за то же время со второго АЦП
                        будет получено kt - 2 отсчетов.*/
        }
        //cout << kt;
        delete noise;
        y3.erase(y3.begin(), y3.end());
        y3.shrink_to_fit();
        y4.erase(y4.begin(), y4.end());
        y4.shrink_to_fit();
        //1.0 ГЕНЕРАЦИЯ МОДЕЛЬНОГО СИГНАЛА
        noise = new double[kt + 1];
        for (int i = 0; i < kt + 1; i++)
        {
            noise[i] = randn();
            //std::cout << noise[i] << endl;

        }
        vector<w_type> yyy;
        for (int i = 0; i < kt + 1; i++)//генерация модельного сигнала (kt+1 вместо kt
        {
            yyy.push_back(sin(2 * M_PI * kp * (i - 1.0) / fi + fi));
            yyy[i] += Q * noise[i];

            //cout << y[i] << endl;
        }
        auto miny = min_element(yyy.begin(), yyy.end(), [](auto a, auto b) {return a.real() < b.real(); });
        auto maxy = max_element(yyy.begin(), yyy.end(), [](auto a, auto b) {return a.real() < b.real(); });
        //cout << endl << miny->real() << ' ' << maxy->real() << endl;
        double key1 = maxy3 - miny3;
        double key2 = maxy4 - miny4;
        ShortComplex* y = new ShortComplex[kt + 1];
        int k = 0;
        for (int i = 0; i < kt + 1; i++)
        {
            y[i].re = yyy[i].real();
            y[i].im = yyy[i].imag();
        }
        yyy.erase(yyy.begin(), yyy.end());
        yyy.shrink_to_fit();
        //2. БПФ разностного сигнала
        universal_fft(y, 199, false); //БПФ

        double* bpf = new double[kt + 1];
        double maxC = 0;
        int kpbpf = 0;
        for (int i = 0; i < kt + 1; i++)
        {
            complex<double> temp = pow(y[i].re, 2) + pow(y[i].im, 2);
            bpf[i] = (temp.real() / (kt + 1.0));
            //cout << bpf[i] << endl;
            //нахождение макс. знач. функции БПФ (целой части количества периодов )
            if (bpf[i] > maxC)
            {
                maxC = bpf[i];
                kpbpf = i;//поиск количества периодов, соответствующих максимуму БПФ
            }
        }
        //cout << kpbpf << " , " << maxC;
        //kpbpf = 10;
        int maxp = 200;
        int maxj = 200;
        int kp20 = kpbpf; //приближенное кол-во периодов разностного сигнала по БПФ
        double* kp2 = new double[maxp];
        double** CC = new double* [maxj];
        double* C = new double[maxp];
        for (int i = 0; i < maxp; i++)
        {
            CC[i] = new double[maxp];
        }
        double* fi_result = new double[maxp];
        int n = 0;
        for (int p = 0; p < maxp; p++)
        {
            kp2[p] = (kpbpf - 1.0) + (2.0 * p) / maxp;
            //cout << kp2[p] << endl;
            double* fi2 = new double[maxj];
            double* skr = new double[maxj];
            double* sskr = new double[maxj];
            for (int j = 0; j < maxj; j++)
            {
                fi2[j] = -M_PI / 2.0 + j * TwoPi / maxj;
                //cout << fi2[j] << endl;
                int kt2 = (kt + 1) * 10; //увеличение кол-ва точек для генерации точного разностного сигнала
                int s = 0;
                int kt2i = kt2 / 10;
                double* z = new double[kt2];
                double* w = new double[kt2i];
                for (int i = 0; i < kt2; i += 10)
                {
                    z[i] = sin(2.0 * M_PI * kp2[p] * (i - 1.0) / kt2 + fi2[j]); //генерация точного разностного сигнала 
                    //cout << z[i] << endl;
                    w[s] = z[i]; //представление разностного сигнала c уменьшенным количеством точек
                    //cout << w[s] << endl;
                    s++;
                }
                delete z;
                skr[j] = 0;
                //cout << "yk and wk" << endl;
                for (int k = 0; k < kt; k++)
                {
                    //cout << y[k].re << " " << w[k] << endl;
                    skr[j] += pow((y[k].re - w[k]), 2);
                }
                /*for (int i = 0; i < kt; i++)
                {
                    cout << w[i] << " " << y[i].re << endl;
                }*/
                delete w;
                sskr[j] = skr[j];
                //cout << sskr[j] << endl;
                CC[j][p] = sskr[j];

            }
            double minSskr = 20000;
            //cout << "sskr = " << endl;
            for (int j = 0; j < maxj; j++)
            {
                //cout << sskr[j] << endl;;

                if (sskr[j] < minSskr)
                {
                    minSskr = sskr[j];
                    n = j; //поиск количества периодов, соответствующих минимуму СКР
                }
            }
            //cout << n << endl;
            //cout << minSskr << "\n\n\n";
            C[p] = minSskr;
            //cout << C[p];
            fi_result[p] = -M_PI / 2.0 + n * TwoPi / maxj; // n - идеально 98.5+
            //cout << fi_result[p] << endl;
            delete fi2;
            delete skr;
            delete sskr;

        }
        int m = 0;
        double minC = 20000;
        for (int i = 0; i < maxp; i++)
        {
            //cout << C[i] << endl;
            if (C[i] < minC)
            {
                m = i;
                minC = C[i];
            }
        }
        //cout << m << "\n\n\n";
        double kp2_res = kpbpf - 1.0 + 2.0 * n / maxp;
        double fi_res = fi_result[m];
        //cout << fi_res;
        if (abs(fi_res) > M_PI)//учет возможного сдвига фазы вторичного сигнала на pi
        {
            fi_res = M_PI - fi_res;
        }
        //cout << fi_res << endl;

        abs_err_grad[j2] = abs(fi - fi_res) / (2.0 * M_PI) * 360.0; //абсолютная погршность определения фазы в градусах

        delete fi_result;
        delete C;
        delete kp2;
        delete CC;
        delete noise;
    } //Цикл для стат.испытаний (конец)
    double aeg_max = -1.0;
    double aeg_std = -1.0;
    double fi_res_std = -1.0;
    double sum = 0.0;
    double mean = 0.0;
    double variance = 0.0;
    //вычисление максимальной и среднеквадратической ошибок
    for (int i = 0; i < maxj2; i++)
    {
        //cout << abs_err_grad[i] << endl;
        if (abs_err_grad[i] > aeg_max)
        {
            aeg_max = abs_err_grad[i];
        }
        sum += abs_err_grad[i];
    }
    mean = sum / maxj2;
    for (int i = 0; i < maxj2; i++)
    {
        variance += pow(abs_err_grad[i] - mean, 2);
    }
    aeg_std = sqrt(variance / maxj2);
    setlocale(LC_ALL, "Russian");
    cout << "Программа завершила " << maxj2 << " стат. испытаний" << endl;
    cout << "Среднеквадратическое отклонение: " << aeg_std << endl;
    cout << "Максимальная ошибка: " << aeg_max << endl;
    delete abs_err_grad;
}