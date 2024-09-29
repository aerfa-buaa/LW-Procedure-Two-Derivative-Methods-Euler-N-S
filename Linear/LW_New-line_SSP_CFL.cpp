#include <math.h>
#include <iostream>
#include <ctime>
#include <cstring>
#include <string>
#include <fstream>
#include <iterator>
#include <valarray>
using namespace std;
// example:line 3-4-L-W   增加CFL数输出误差   Line Excample-1  精简版本只有1up,2cneter   本程序方程U_t-U_x=0   验证了与U_t+U_x=0一致
// g++ -O2 -std=c++11 -Wall "-Wl,--stack=268435456" TDMS-line_SSP_CFL.cpp -o run.exe

// g++ -O2 -std=c++11 -o run.exe TDMS-line_SSP_CFL.cpp
class cfda
{

private:
    int w;

public:
    double PI = 3.141592653589793;
    // double PI = 3.1415926535897932385;

    int n;                                         // 网格结点数（有限差分）
    int nb;                                        // 边界网格结点数（有限差分）
    double t, tt, CFL;                             // 计算总时 间t,实际计算时间tt;
    int nt, kt, ij, ijj, onet_out = 2;             // 计算时间步长;
    double nx_begin, nx_end;                       // 网格长度;
    int n1;                                        // 总的网格结点数（有限差分）
    double *Ltt_n0, *x;                            // 声明变长数组
    double *L_n0, *L_n1, *L_n2, *L_n3, *L_n4;      // f=L
    double *Lt_n0, *Lt_n1, *Lt_n2, *Lt_n3, *Lt_n4; // L_t

    double *TL_n0, *TL_n1, *TL_n2, *TL_n3, *TL_n4;      // Tf=L
    double *TLt_n0, *TLt_n1, *TLt_n2, *TLt_n3, *TLt_n4; // TL_t

    double *u_n0, *u_n1, *u_n2, *u_n3, *u_n4, *u_n5, *u_n6, *u_n7, *u_n8, *u_n9, *f_k; // 推进步中间时刻
    double *f_n0, *f_n1, *f_n2, *f_n3, *f_n4, *f_n5, *f_n6, *f_n7, *f_n8, *f_n9;       // 推进步中间时刻
    double *u_nn, *f_nn, *Lttx_n0;                                                     // 下一时刻（n+1)
    double *Tu1, *Tu2;
    // double  u_excat[513 + 6 * 2];
    double *df, *u_excat; // du[513 + 6 * 2],df[513 + 6 * 2],
    ////-----------------------------------------
    double dx;
    double dt;
    double *dd; // ss
    double A1, A2, A3, A4, A5, b0, b1, b2, b3, B0, B1, B2, bb1, bb2, bbb1, bbb2, C0, C1, C2, C3, D0, D1, D2, D3, maxuu, TVmax, k;
    double a21, a32, a31, aa12, aa21, aa31, aa32, w1, w2, v0, v1, v2, v3, vv0, vv1, vv2, vv3, ww1, ww2, d31, d32, a22, aa22;
    double put_obj[2][1300];   //
    double put_one_n[2][1300]; // onet_out = 10;
    ////-----------------------------------------
    void intc()
    {
        dx = (nx_end - nx_begin) / (1.0 * (n - 1));

        // for (int i = 0; i < n1; i++)
        // {
        //     x[i] = nx_begin + (i - nb) * dx;
        //     u_n0[i] = intc_fun(x[i]);
        //     df[i] = 1;
        // }

        for (int i = 0; i < n1; i++) // 间断
        {
            x[i] = nx_begin + (i - nb) * dx;
            u_n0[i] = 0.0;
            if (x[i] >= -0.0 && x[i] <= 0.5)
            {
                u_n0[i] = 1.0;
            }
            // df[i] = 1;
        }
        // Write_obj(0);
    }

    double intc_fun(double xxf)
    {
        double yyf;
        // yyf = (sin(2 * PI * xxf) / 2.0 + 0.5) * 1.0; //  / PIbergure
        // yyf = sin(PI * xxf) / 2.0 + 0.25;//bergure
        yyf = sin(PI * xxf); // line
        return yyf;
    }

    // //-----------------------------------------
    void border()
    {
        borderfun(u_n0), borderfun(u_n1), borderfun(u_n2), borderfun(u_n3), borderfun(u_n4), borderfun(u_n5);
        borderfun(u_n6), borderfun(u_n7), borderfun(u_n8), borderfun(u_n9), borderfun(u_nn);
        borderfun(f_n0), borderfun(f_n1), borderfun(f_n2), borderfun(f_n3), borderfun(f_n4), borderfun(f_n5);
        borderfun(f_n6), borderfun(f_n7), borderfun(f_n8), borderfun(f_n9), borderfun(f_nn);
        borderfun(L_n0), borderfun(L_n1), borderfun(L_n2), borderfun(L_n3), borderfun(L_n4);
        borderfun(Lt_n0), borderfun(Lt_n1), borderfun(Lt_n2), borderfun(Lt_n3), borderfun(Lt_n4);
        borderfun(Ltt_n0), borderfun(Lttx_n0);
    }

    void borderfun(double *ffbc)
    {

        for (int i = 0; i < nb; i++)
        {

            *(ffbc + i) = *(ffbc + n1 - 2 * nb + i - 1);  // 周期边条
            *(ffbc + n1 - nb + i) = *(ffbc + nb + i + 1); // 周期边条

            // *(ffbc + i) = 0.;
            // *(ffbc + n1 - nb + i) = 0.;

            // f_n0[i] = f_n0[nb];                   //恒定边条
            // f_n0[n1 - 1 - i] = f_n0[n1 - 1 - nb]; //恒定边条

            // *(ffbc + i) = *(ffbc + nb);                   //恒定边条
            // *(ffbc + n1 - 1 - i) = *(ffbc + n1 - 1 - nb); //恒定边条

            // *(ffbc + i) = 0.0;                   //恒定边条
            // *(ffbc + n1 - 1 - i) = 0.0; //恒定边条
            // ffbc++;
        }
    }
    // //-----------------------------------------
    void carr(double *p1, double *p2, int i)
    {
        while (i-- > 0)
        {
            *p1++ = *p2++;
        }
    }

    // //-----------------------------------------
    void compt_2_stage_int()
    {
        int i, j;
        a21 = 0.765; // 0.532  0.468; //0.765;   // 1. / 210. * (93. + sqrt(2139.));
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);

            TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void TDRK23_line_SSP()
    {

        int i, j;

        // a21 = 0.594223212099088, v2 = 0.306027487008159;
        a21 = 2. / 3., v2 = 3. / 8.; // Tylor
        ww1 = 0.0, ww2 = 0.0, w1 = 0.0, w2 = 0.0;
        v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            // Lt_n0[j] = -1. * dfdxcent(L_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + w1 * TL_n0[j] + w2 * TL_n1[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + ww1 * TLt_n0[j] + ww2 * TLt_n1[j]);

            // TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            // TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void TDRK34_line_SSP()
    {
        int i, j;
        double uu, vv;
        a21 = 0.532230958301739;
        a31 = 0.457385979790557;
        a32 = 0.207902718086617;
        aa31 = 0.0553261314403881;
        v3 = 0.374215233358609;

        aa21 = pow(a21, 2) / 2.;
        aa32 = (pow(a31, 2) - 2. * a21 * a32 + 2. * a31 * a32 + pow(a32, 2) - 2. * aa31) / 2.;
        uu = a31 + a32;
        v2 = -(((3. * pow(a21, 2) * a32 + 6. * a21 * aa31 - 3. * a21 * pow(uu, 2) + pow(uu, 3)) * v3) / pow(a21, 3));
        v1 = 1. - v2 - v3;
        vv1 = -(-1. + 2. * uu - 2. * a21 * (-1. + 3. * uu + a21 * (a21 - 3. * uu) * v2) + 2. * (3. * a21 - uu) * pow(uu, 2) * v3) / (12. * a21 * uu);
        vv2 = -(-1. + 2. * uu + 4. * pow(a21, 3) * v2 - 6. * pow(a21, 2) * uu * v2 - 2. * pow(uu, 3) * v3) / (12. * a21 * (a21 - uu));
        vv3 = -(-1. + 2. * a21 - 2. * pow(a21, 3) * v2 - 6. * a21 * pow(uu, 2) * v3 + 4. * pow(uu, 3) * v3) / (12. * uu * (-a21 + uu));

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = u_n0[j] + dt * (a31 * L_n0[j] + a32 * L_n1[j]) + dt * dt * (aa31 * Lt_n0[j] + aa32 * Lt_n1[j]);
            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = dfdx(f_n2, j);
            Lt_n2[j] = dfdxx(f_n2, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + v3 * L_n2[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + vv3 * Lt_n2[j]);

            // TL_n1a[j] = L_n1[j], TL_n1[j] = L_n0[j];
            // TLt_n1a[j] = Lt_n1[j], TLt_n1[j] = Lt_n0[j];
        }
    }

    // //-----------------------------------------

    void compt_TDMS_stage_int()
    {
        int i, j;
        i = n1 - nb * 2;
        // border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);

            Tu2[j] = Tu1[j], TL_n1[j] = TL_n0[j], TLt_n1[j] = TLt_n0[j];
            Tu1[j] = u_n0[j], TL_n0[j] = L_n0[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt_TDMS34_line()
    {
        int i, j; // ynn= -b0 *yn + dt* v0 *L + dt^2 * vv0*Lt;    Tu1=u^(n-1);  Tu2=u^(n-2);
        double uu, vv;
        b1 = 0.284944540740210;
        b2 = 0.156281668047679;
        v1 = 0.491862931463468;

        w1 = 1.;
        w2 = 1.;
        b3 = 1. - b1 - b2;
        vv2 = 0.;
        uu = w1 + w2, vv = (b1 + b2) * w2;
        vv1 = (3. + 12. * pow(uu, 2) * v1 * pow(w1, 2) - b1 * pow(uu, 3) * (3. * w1 - w2) + 4. * w2 + b2 * pow(w2, 4)) / (12. * uu * pow(w1, 2) * (3. * w1 + w2));

        v2 = (1. + w2 + w1 * (1. + pow(uu, 3) * v1 - 6. * pow(uu, 2) * vv1 * w1 + b2 * pow(w2, 3))) / ((3. * w1 - w2) * pow(w2, 3));

        vv3 = (3. + (b1 - 12. * vv1) * pow(w1, 4) + 2. * (b1 - 6. * vv1) * pow(w1, 3) * w2 + 2. * w1 * (2. + 3. * w2 - vv * pow(w2, 2)) + w2 * (8. + 6. * w2 - vv * pow(w2, 2))) / (12. * uu * w2);

        v3 = (-1. - 2. * uu - (b1 - 2. * v1) * pow(w1, 4) + 2. * (3. + 2. * vv) * w1 * pow(w2, 2) + (2. + vv) * pow(w2, 3) + pow(w1, 3) * (-4. * b1 * w2 + 6. * v1 * w2)) / (2. * pow(w2, 2) * (3. * w1 + w2));

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_nn[j] = b1 * Tu2[j] + b2 * Tu1[j] + b3 * u_n0[j] +
                      dt * (v1 * TL_n1[j] + v2 * TL_n0[j] + v3 * L_n0[j]) +
                      dt * dt * (vv1 * TLt_n1[j] + vv2 * TLt_n0[j] + vv3 * Lt_n0[j]);

            Tu2[j] = Tu1[j], TL_n1[j] = TL_n0[j], TLt_n1[j] = TLt_n0[j];
            Tu1[j] = u_n0[j], TL_n0[j] = L_n0[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt_TDMS23_line()
    {
        int i, j; // ynn= -b0 *yn + dt* v0 *L + dt^2 * vv0*Lt;    Tu1=u^(n-1);  Tu2=u^(n-2);

        // b1 = 0.2;
        b1 = 0.416407864998738;
        w1 = 1.0;

        b2 = 1.0 - b1;
        vv1 = 0.0;
        v1 = 1. / 3. * (b1 + 1. / pow(w1, 3));
        v2 = 1. - 1. / (3.0 * w1 * w1) + (2.0 * b1 * w1) / 3.0;
        vv2 = (2. + 3.0 * w1 - b1 * pow(w1, 3)) / (6.0 * w1);

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_nn[j] = b1 * Tu1[j] + b2 * u_n0[j] + dt * (v1 * TL_n0[j] + v2 * L_n0[j]) +
                      dt * dt * (vv1 * TLt_n0[j] + vv2 * Lt_n0[j]);

            Tu1[j] = u_n0[j], TL_n0[j] = L_n0[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    ////-----------------------------------------

    void RK3_compt_t1()
    {
        // compt_t1(cff1.uu_t1, cff1.uu_t0,n_all);
        // uu_t1 = uu_t0 + dt * df;  df_dx(double *g, int i)
        int j = 0, i;
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            // df[nb + j] = -df_dx(f_n0, nb + j); //-uu_t0[nb + j]; -dfdx
            df[nb + j] = dfdx(f_n0, nb + j);
            u_n1[nb + j] = u_n0[nb + j] + dt * df[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        j = 0, i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            // df[nb + j] = -df_dx(f_n1, nb + j); //-uu_t1[nb + j];
            df[nb + j] = dfdx(f_n1, nb + j);
            u_n2[nb + j] = u_n0[nb + j] * 0.75 + u_n1[nb + j] * 0.25 + 0.25 * dt * df[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        j = 0, i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            // df[nb + j] = -df_dx(f_n2, nb + j);
            df[nb + j] = dfdx(f_n2, nb + j);
            u_nn[nb + j] = u_n0[nb + j] * 1.0 / 3.0 + u_n2[nb + j] * 2.0 / 3.0 + 2.0 / 3.0 * dt * df[nb + j];
            j++;
        }
    }

    void SSPRK54_compt()
    {
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 0.39175222700392 * dt * df[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = 0.44437049406734 * u_n0[nb + j] + 0.55562950593266 * u_n1[nb + j] + 0.36841059262959 * dt * df[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = dfdx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = 0.62010185138540 * u_n0[nb + j] + 0.37989814861460 * u_n2[nb + j] + 0.25189177424738 * dt * df[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = dfdx(f_n3, nb + j); //-uu_t1[nb + j];
            Ltt_n0[nb + j] = df[nb + j];
            u_n4[nb + j] = 0.17807995410773 * u_n0[nb + j] + 0.82192004589227 * u_n3[nb + j] + 0.54497475021237 * dt * df[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }

        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = dfdx(f_n4, nb + j);
            u_nn[nb + j] = 0.00683325884039 * u_n0[nb + j] + 0.51723167208978 * u_n2[nb + j] + 0.12759831133288 * u_n3[nb + j] +
                           0.34833675773694 * u_n4[nb + j] + 0.08460416338212 * dt * Ltt_n0[nb + j] + 0.22600748319395 * dt * df[nb + j];
            j++;
        }
    }

    //------------------------------------------------
    void f_eq_u(double *ff, double *uu)
    {
        for (int i = 0; i < n1; i++)
        {

            *ff++ = f_u(*uu++);
        }
    }
    double f_u(double xx)
    {
        double y;
        // y = xx * xx / 2.0;//bergure
        y = xx; // line
        return y;
    }
    void Store_obj(double *ss, double *ff)
    {
        for (int ik = 0; ik < n1; ik++) //
        {
            *ss++ = *ff++;
        }
    }

    void Write_obj(int j)
    {
        if (j < -0.9)
        {
            // string title = "result_t_burgers_excat"; //+ to_string(t);
            string title = to_string(n);
            string Title = title + ".plt";
            ofstream ofs;
            ofs.open(Title, ios::out);
            ofs << " VARIABLES=x,u,hy,hy0,fhy " << endl;

            for (int i = nb; i < n1 - nb; i++)
            {

                ofs.precision(10);
                ofs << x[i];
                for (int j = 0; j < kt; j++)
                {
                    ofs << "  ";
                    ofs.precision(18);
                    ofs << put_obj[j][i];
                }
                ofs << "  " << f_k[i] << "  " << f_n7[i] << "  " << f_k[i] << endl;
            }
            ofs.close();
        }

        if (j < 0.5 && j > -0.8)
        {
            // string title = "result_t_burgers_excat"; //+ to_string(t);
            string title = to_string(n);
            string Title = "excat.txt";
            ofstream ofs;
            ofs.open(Title, ios::out);
            ofs << "Title=" << title << endl;

            for (int i = nb; i < n1 - nb; i++)
            {

                ofs.precision(10);
                ofs << x[i];
                ofs << "  ";
                ofs.precision(18);
                ofs << *(u_excat + i);
                ofs << endl;
            }
            ofs.close();
        }

        if (j > 1.5)
        {
            string title = "result_one_t_burgers_every-" + to_string(j);
            string Title = title + ".plt";
            ofstream ofs;
            ofs.open(Title, ios::out);
            ofs << "Title=" << title << endl;
            for (int i = nb; i < n1 - nb; i++)
            {

                ofs.precision(10);
                ofs << x[i];
                for (j = 0; j < onet_out; j++)
                {
                    ofs << "  ";
                    ofs.precision(18);
                    ofs << put_one_n[j][i];
                }

                ofs << endl;
            }
            ofs.close();
        }
    }

    //---------------------------------------------------------

    int Perimeter()
    {
        cout << " yes ";
        return 2 * t;
    }

    double dfdx(double *g, int i)
    {
        double y;
        // 1up
        // y = (*(g + i) - *(g + i - 1)) / dx; //u_t+u_x=0
        y = (*(g + i + 1) - *(g + i)) / dx; // u_t-u_x=0
        return y;
    }

    double dfdxx(double *g, int i)
    {
        double y;
        // 2ord
        // y = (*(g + i + 1) - *(g + i) * 2.0 + *(g + i - 1)) / dx / dx; // u_t-u_x=0
        y = (*(g + i + 2) - *(g + i + 1) * 2.0 + *(g + i)) / dx / dx; // u_t-u_x=0
        // y = (*(g + i ) - *(g + i - 1) * 2.0 + *(g + i-2)) / dx / dx; //u_t+u_x=0
        return y;
    }

 

    void compt_CFL()
    {
        // double *umax;
        // umax = max_element(u_n0 + 1, u_n0 + n1 - 1);
        // cout << *umax;
        // dt = CFL * dx / (*umax);
        // dt=0.02/1.0;
        dt = CFL * dx;
    }

    void maxTV()
    {
        double TV = 0.0;
        for (int i = 0; i < n1 - 1; i++) // 间断
        {
            TV = abs(u_nn[i + 1] - u_nn[i]) + TV;
        }
        TV = abs(TV - 2.0);
        // cout << " \n TV: " << TV;
        if (tt > 0.0)
        {
            if (TV > TVmax)
            {
                TVmax = TV;
            }
        }
    }

    void TVabs(int m)
    {
        double TV = 0.0;
        TV = TVmax;

        cout << "\n  CFL: " << CFL << ",  TV: " << TV << "   \n";
        string Title = "TV_CFL.plt";
        ofstream ofs;
        if (m < 2)
        {
            ofs.open(Title, ios::out);
            ofs << " VARIABLES=CFL,TV " << endl;
            ofs.precision(6), ofs << CFL, ofs << "  ";
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
        else
        {
            ofs.open(Title, ios::app);
            ofs.precision(6), ofs << CFL, ofs << "  ";
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
    }

    void TVabs_t(int m)
    {
        double TV = 0.0;
        TV = TVmax;

        if (TV < 2.0E-16)
        {
            TV = 2.0E-16;
        }
        // cout << "\n  CFL: " << CFL << ",  TV: " << TV << "   \n";
        string Title = "TV_t.plt";
        ofstream ofs;
        if (m < 2)
        {
            ofs.open(Title, ios::out);
            ofs << " VARIABLES=t,TV " << endl;
            ofs << m, ofs << "  ";
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
        else
        {
            ofs.open(Title, ios::app);
            ofs << m, ofs << "  ";
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
    }

    void compt_shock()
    {
        int i, j;
        double date0, eps, df1, df2, ka;
        date0 = 0.6;
        eps = 0.9 * date0 / (1. - 0.9 * date0);
        i = n1 - nb * 2;
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            df1 = f_n0[j + 1] - f_n0[j], df2 = f_n0[j] - f_n0[j - 1];
            // f_k[j] = (2. * abs(df1 * df2) + eps) / (df1 * df1 + df2 * df2 + eps);  )
            f_n7[j] = abs(df1) + abs(df2);
            f_k[j] = 0.;
        }
        borderfun(f_k);

        i = n1 - nb * 2;
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            for (int k = -2; k < 3; k++)
                if (f_n7[j + k] > 1E-14)
                    f_k[j] = 1.;
        }
        borderfun(f_k);
    }

    void Write_hy(int m)
    {

        string Title = "hy.plt";
        ofstream ofsy;
        if (m < 1)
        {
            ofsy.open(Title, ios::out);
            ofsy << " VARIABLES=x,t " << endl;
            for (int i = nb; i < n1 - nb; i++) // 间断  int i = nb; i < n1 - nb; i++
            {
                if (f_k[i] > 0.01)
                    ofsy << x[i] << "  " << m << endl;
            }
            ofsy.close();
        }
        else
        {
            ofsy.open(Title, ios::app);
            for (int i = nb; i < n1 - nb; i++) // 间断  int i = nb; i < n1 - nb; i++
            {
                if (f_k[i] > 0.01)
                    ofsy << x[i] << "  " << m << endl;
            }
            ofsy.close();
        }
    }
};
////-----------------------------------------
void comput_tt(cfda &);
int comput_main(int, double, double, int);
int main()
{
    // int n, oder_n = 2, begin_cell = 1600;
    int n, oder_n = 139, begin_cell = 3200;
    double t = 100.0, CFL = 0.01; // 计算时间总长，CFL数  CFL = 0.66 1.45
    n = begin_cell + 1;
    for (int m = 1; m < oder_n; m++)
    {
        comput_main(n, t, CFL, m);
        CFL = CFL + 0.02;
    }
    // system("pause");
    return 2;
}
int comput_main(int n, double t, double CFL, int m)
{
    int nb = 6; // 边界网格结点数（有限差分）
    int n1 = n + 6 * 2;
    double main_df[n + 6 * 2], u_excat[n + 6 * 2];
    double Ltt_n0[n + 6 * 2], x[n + 6 * 2];                                                          // 声明变长数组
    double L_n0[n + 6 * 2], L_n1[n + 6 * 2], L_n2[n + 6 * 2], L_n3[n + 6 * 2], L_n4[n + 6 * 2];      // f=L
    double Lt_n0[n + 6 * 2], Lt_n1[n + 6 * 2], Lt_n2[n + 6 * 2], Lt_n3[n + 6 * 2], Lt_n4[n + 6 * 2]; // L_t
    double u_n0[n + 6 * 2], u_n1[n + 6 * 2], u_n2[n + 6 * 2], u_n3[n + 6 * 2], u_n4[n + 6 * 2];      // 推进步中间时刻
    double u_n5[n + 6 * 2], u_n6[n + 6 * 2], u_n7[n + 6 * 2], u_n8[n + 6 * 2], u_n9[n + 6 * 2];      // 推进步中间时刻
    double f_n0[n + 6 * 2], f_n1[n + 6 * 2], f_n2[n + 6 * 2], f_n3[n + 6 * 2], f_n4[n + 6 * 2];      // 推进步中间时刻
    double f_n5[n + 6 * 2], f_n6[n + 6 * 2], f_n7[n + 6 * 2], f_n8[n + 6 * 2], f_n9[n + 6 * 2];      // 推进步中间时刻
    double u_nn[n + 6 * 2], f_nn[n + 6 * 2], Lttx_n0[n + 6 * 2], f_k[n + 6 * 2];                     // 下一时刻（n+1)

    double TL_n0[n + 6 * 2], TL_n1[n + 6 * 2], TL_n2[n + 6 * 2], TL_n3[n + 6 * 2], TL_n4[n + 6 * 2];
    double TLt_n0[n + 6 * 2], TLt_n1[n + 6 * 2], TLt_n2[n + 6 * 2], TLt_n3[n + 6 * 2], TLt_n4[n + 6 * 2];
    double Tu1[n + 6 * 2], Tu2[n + 6 * 2];
    class cfda cff;
    cff.kt = 1;
    cff.A1 = 1 / 3.0;
    cff.A2 = 2 / 3.0;
    cff.t = t;           // 1.0 / cff.PI; //计算总时间
    cff.nx_begin = -1.0; // 网格右端点
    cff.nx_end = 1.0;    // 网格左端点
    cff.CFL = CFL, cff.TVmax = 0.0;

    cff.n = n, cff.nb = nb, cff.n1 = n1;
    cff.Ltt_n0 = Ltt_n0, cff.x = x;
    cff.u_n0 = u_n0, cff.u_n1 = u_n1, cff.u_n2 = u_n2, cff.u_n3 = u_n3, cff.u_n4 = u_n4;                     // 推进步中间时刻
    cff.u_n5 = u_n5, cff.u_n6 = u_n6, cff.u_n7 = u_n7, cff.u_n8 = u_n8, cff.u_n9 = u_n9;                     // 推进步中间时刻
    cff.f_n0 = f_n0, cff.f_n1 = f_n1, cff.f_n2 = f_n2, cff.f_n3 = f_n3, cff.f_n4 = f_n4;                     // 推进步中间时刻
    cff.f_n5 = f_n5, cff.f_n6 = f_n6, cff.f_n7 = f_n7, cff.f_n8 = f_n8, cff.f_n9 = f_n9, cff.f_k = f_k;      // 推进步中间时刻
                                                                                                             // 声明变长数组
    cff.u_nn = u_nn, cff.f_nn = f_nn, cff.Lttx_n0 = Lttx_n0;                                                 // 下一时刻（n+1)
    cff.L_n0 = L_n0, cff.L_n1 = L_n1, cff.L_n2 = L_n2, cff.L_n3 = L_n3, cff.L_n4 = L_n4;                     // f=L
    cff.Lt_n0 = Lt_n0, cff.Lt_n1 = Lt_n1, cff.Lt_n2 = Lt_n2, cff.Lt_n3 = Lt_n3, cff.Lt_n4 = Lt_n4;           // L_t
                                                                                                             // 存储一时刻（n-1)
    cff.TL_n0 = TL_n0, cff.TL_n1 = TL_n1, cff.TL_n2 = TL_n2, cff.TL_n3 = TL_n3, cff.TL_n4 = TL_n4;           // f=TL
    cff.TLt_n0 = TLt_n0, cff.TLt_n1 = TLt_n1, cff.TLt_n2 = TLt_n2, cff.TLt_n3 = TLt_n3, cff.TLt_n4 = TLt_n4; // TL_t

    cff.Tu1 = Tu1, cff.Tu2 = Tu2;
    cff.df = main_df, cff.u_excat = u_excat, cff.intc();
    for (cff.ij = 0; cff.ij < cff.kt; cff.ij++)
    {
        // cff.nt = 1000 * pow(2, cff.ij); //计算时间步数

        comput_tt(cff);

        cff.Store_obj(cff.put_obj[cff.ij], cff.u_nn); // 存储数据
        // cff.Write_obj(cff.ij + 1);                    //存储数据
    }
    cff.Write_obj(-1);
    cff.TVabs(m);
    return 2;
}
//-----------------------------------------
void comput_tt(cfda &cff1)
{
    int n_all, i3 = 0, tt_flag = 0;
    cff1.tt = 0, n_all = cff1.n1, cff1.border();
    cff1.f_eq_u(cff1.f_n0, cff1.u_n0), cff1.TVmax = 0.0;
    for (int i1 = 0; i1 <= 100; i1++) // 500 200 50
    {
        cff1.compt_CFL();
        cff1.tt = cff1.tt + cff1.dt;
        if (i1 < 4) // 4  5
        {
            cff1.compt_TDMS_stage_int();

            cff1.TDRK23_line_SSP();
            // cff1.TDRK34_line_SSP();

            // cff1.RK3_compt_t1(); // 3-RK
            // cff1.SSPRK54_compt(); //
        }
        else
        {
            cff1.compt_TDMS23_line();
            // cff1.compt_TDMS34_line();
        }
        //___________________________________________________________________
        if (cff1.tt > cff1.t)
        {
            cff1.f_eq_u(cff1.f_nn, cff1.u_nn);
            cout.precision(18); // 精度为18，正常为6
            cout << " \n comput time: " << cff1.tt - cff1.dt << " comput time n:   " << i1;
            break;
        }
        cff1.carr(cff1.u_n0, cff1.u_nn, n_all);
        cff1.f_eq_u(cff1.f_n0, cff1.u_nn);

        // cout << " \n comput time: " << cff1.tt << " comput time n:   " << i1;

        cff1.border();
        cff1.maxTV();
        cff1.TVabs_t(i1);

        if (i3 * cff1.t / cff1.onet_out < cff1.tt + 0.00000001)
        {
            // cff1.Store_obj(cff1.put_one_n[i3], cff1.u_nn); //存储数据
            i3++;
        }
    }
}
