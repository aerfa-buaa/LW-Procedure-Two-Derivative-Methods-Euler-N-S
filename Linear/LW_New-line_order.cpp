#include <math.h>
#include <iostream>
#include <ctime>
#include <cstring>
#include <string>
#include <fstream>
#include <iterator>
#include <valarray>
using namespace std;
// example:line 3-4-L-W   测试精度 算例3 weno, 不计算Tv  简化函数  计算误差值  U_t-U_x=0
class cfda
{

private:
    int w;

public:
    double PI = 3.141592653589793;
    // double PI = 3.1415926535897932385;

    int n;                                         // 网格结点数（有限差分）
    int nb;                                        // 边界网格结点数（有限差分）
    double t, tt, CFL, xt, xt1;                    // 计算总时 间t,实际计算时间tt;
    int nt, kt, ij, ijj, onet_out = 2;             // 计算时间步长;
    double nx_begin, nx_end;                       // 网格长度;
    int n1;                                        // 总的网格结点数（有限差分）
    double *Ltt_n0, *x;                            // 声明变长数组
    double *L_n0, *L_n1, *L_n2, *L_n3, *L_n4;      // f=L
    double *Lt_n0, *Lt_n1, *Lt_n2, *Lt_n3, *Lt_n4; // L_t

    double *TL_n0, *TL_n1, *TL_n2, *TL_n3, *TL_n4;      // Tf=L
    double *TLt_n0, *TLt_n1, *TLt_n2, *TLt_n3, *TLt_n4; // TL_t

    double *u_n0, *u_n1, *u_n2, *u_n3, *u_n4, *u_n5, *u_n6, *u_n7, *u_n8, *u_n9; // 推进步中间时刻
    double *f_n0, *f_n1, *f_n2, *f_n3, *f_n4, *f_n5, *f_n6, *f_n7, *f_n8, *f_k;  // 推进步中间时刻
    double *u_nn, *f_nn, *Lttx_n0;                                               // 下一时刻（n+1)
    double *Tu1, *Tu2;
    // double  u_excat[513 + 6 * 2];
    double *df, *u_excat; // du[513 + 6 * 2],df[513 + 6 * 2],
    ////-----------------------------------------
    double dx;
    double dt;
    double *dd; // ss
    double A1, A2, A3, AA1, AA2, AA3, b0, b1, b2, b3, B0, B1, B2, B3, bb1, bb2, bb3, bbb1, bbb2, bbb3, C0, C1, C2, C3, D0, D1, D2, D3, maxuu, TVmax, k;
    double a21, a22, a31, a32, aa12, aa21, aa22, aa31, aa32, aaa31, aaa32, w1, d31, d32,
        w2, v0, v1, v2, v3, vv0, vv1, vv2, vv3, ww1, ww2;
    double put_obj[2][13000];   //
    double put_one_n[2][13000]; // onet_out = 10;

    ////-----------------------------------------
    void intc()
    {
        dx = (nx_end - nx_begin) / (1.0 * (n - 1.0));
        for (int i = 0; i < n1; i++)
        {
            x[i] = nx_begin + (i - nb) * dx;
            u_n0[i] = intc_fun(x[i]);
            df[i] = 1;
        }

        // for (int i = 0; i < n1; i++) //间断
        // {
        //     x[i] = nx_begin + (i - nb) * dx;
        //     u_n0[i] = 0.1;
        //     if (x[i] >= 0.1 && x[i] <= 0.5)
        //     {
        //         u_n0[i] = 1.0;
        //     }
        //     // df[i] = 1;
        // }
        // Write_obj(0);
    }

    double intc_fun(double xxf)
    {
        double yyf;
        // yyf = (sin(2 * PI * xxf) / 2.0 + 0.5) * 1.0; //  / PI bergure
        // yyf = sin(PI * xxf) / 2.0 + 0.25;//bergure
        yyf = 0.5 * sin(PI * xxf) + 0.5; // line
        // yyf = 0.5 * sin(PI * xxf) + 0.5; //line
        return yyf;
    }

    // //-----------------------------------------
    void border()
    {
        borderfun(u_n0), borderfun(u_n1), borderfun(u_n2), borderfun(u_n3), borderfun(u_n4), borderfun(u_n5);
        borderfun(u_n6), borderfun(u_n7), borderfun(u_n8), borderfun(u_n9), borderfun(u_nn);
        borderfun(f_n0), borderfun(f_n1), borderfun(f_n2), borderfun(f_n3), borderfun(f_n4), borderfun(f_n5);
        borderfun(f_n6), borderfun(f_n7), borderfun(f_n8), borderfun(f_k), borderfun(f_nn);
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

            // *(ffbc + i) = 0;
            //*(ffbc + n1 - nb + i) = 0;
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

    //------------------------------------------------
    void compt23_LLttnn_line()
    {
        int i, j;
        // a21 = 0.594223212099088, v2 = 0.306027487008159;
        a21 = 2. / 3., v2 = 3. / 8.;
        ww1 = 0.0, ww2 = 0.0, w1 = 0.0, w2 = 0.0;
        v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = df_dx_up(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            // Lt_n0[j] = dfdxx(f_n0, j);
            Lt_n0[j] = dfdxcent(L_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = df_dx_up(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            // Lt_n1[j] = dfdxx(f_n1, j);
            Lt_n1[j] = dfdxcent(L_n1, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + w1 * TL_n0[j] + w2 * TL_n1[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + ww1 * TLt_n0[j] + ww2 * TLt_n1[j]);

            // TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            // TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt34_LLttnn_line()
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
            L_n0[j] = df_dx_up(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxcent(L_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = df_dx_up(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxcent(L_n1, j);
            u_n2[j] = u_n0[j] + dt * (a31 * L_n0[j] + a32 * L_n1[j]) + dt * dt * (aa31 * Lt_n0[j] + aa32 * Lt_n1[j]);
            f_n2[j] = f_u(u_n2[j]);
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = df_dx_up(f_n2, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n2[j] = dfdxcent(L_n2, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + v3 * L_n2[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + vv3 * Lt_n2[j]);

            // TL_n1a[j] = L_n1[j], TL_n1[j] = L_n0[j];
            // TLt_n1a[j] = Lt_n1[j], TLt_n1[j] = Lt_n0[j];
        }
    }

    //----------------------------------------
    void compt_TDMS_23_ALW_VSS()
    {
        int i, j;
        w1 = xt / dt;
        // b1 = 0.2;
        // b1 = 0.2013 * pow(w1, -1.592);
        b1 = 0.4091 * pow(w1, -1.861); // h_b1 = 0.;  //Tylor
        b2 = 1.0 - b1;
        vv1 = 0.0;
        v1 = 1. / 3. * (b1 + 1. / pow(w1, 3));
        v2 = 1. - 1. / (3.0 * w1 * w1) + (2.0 * b1 * w1) / 3.0;
        vv2 = (2. + 3.0 * w1 - b1 * pow(w1, 3)) / (6.0 * w1);

        border();
        i = n1 - nb * 2;
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = df_dx_up(f_n0, j);
        }
        border();
        i = n1 - nb * 2;
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxcent(L_n0, j);
            u_nn[j] = b1 * Tu1[j] + b2 * u_n0[j] + dt * (v1 * w1 * TL_n0[j] + v2 * L_n0[j]) +
                      dt * dt * (vv1 * w1 * w1 * TLt_n0[j] + vv2 * Lt_n0[j]);

            Tu1[j] = u_n0[j], TL_n0[j] = L_n0[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt_TDMS_34_ALW_VSS()
    {
        int i, j;
        double uu, vv;
        w2 = xt / dt, w1 = xt1 / dt;

        b1 = -0.283 + 0.088 * w1 * w1 + 6.246 * exp(-1.227 * w1 - 1.354 * w2);
        b2 = 0.205 + 0.3156 * pow(w1, -1.508) * pow(w2, 1.065) - 3.1 * exp(-2.123 * w1);

        v1 = 22.23 + 8.837 * w1 * w1 - 14.05 * exp(0.7255 * w1 - 0.1072 * w2) - 9.275 * sin(0.5058 * w2);

        b3 = 1. - b1 - b2;
        vv2 = 0.;
        uu = w1 + w2, vv = (b1 + b2) * w2;
        vv1 = (3. + 12. * pow(uu, 2) * v1 * pow(w1, 2) - b1 * pow(uu, 3) * (3. * w1 - w2) + 4. * w2 + b2 * pow(w2, 4)) / (12. * uu * pow(w1, 2) * (3. * w1 + w2));

        v2 = (1. + w2 + w1 * (1. + pow(uu, 3) * v1 - 6. * pow(uu, 2) * vv1 * w1 + b2 * pow(w2, 3))) / ((3. * w1 - w2) * pow(w2, 3));

        vv3 = (3. + (b1 - 12. * vv1) * pow(w1, 4) + 2. * (b1 - 6. * vv1) * pow(w1, 3) * w2 + 2. * w1 * (2. + 3. * w2 - vv * pow(w2, 2)) + w2 * (8. + 6. * w2 - vv * pow(w2, 2))) / (12. * uu * w2);

        v3 = (-1. - 2. * uu - (b1 - 2. * v1) * pow(w1, 4) + 2. * (3. + 2. * vv) * w1 * pow(w2, 2) + (2. + vv) * pow(w2, 3) + pow(w1, 3) * (-4. * b1 * w2 + 6. * v1 * w2)) / (2. * pow(w2, 2) * (3. * w1 + w2));

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = df_dx_up(f_n0, j);
        }

        border();
        i = n1 - nb * 2;
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxcent(L_n0, j);
            u_nn[j] = b1 * Tu2[j] + b2 * Tu1[j] + b3 * u_n0[j] +
                      dt * (v1 * w1 * TL_n1[j] + v2 * w2 * TL_n0[j] + v3 * L_n0[j]) +
                      dt * dt * (vv1 * w1 * w1 * TLt_n1[j] + vv2 * w2 * w2 * TLt_n0[j] + vv3 * Lt_n0[j]);

            Tu2[j] = Tu1[j], TL_n1[j] = TL_n0[j], TLt_n1[j] = TLt_n0[j];
            Tu1[j] = u_n0[j], TL_n0[j] = L_n0[j], TLt_n0[j] = Lt_n0[j];
        }
    }

     ////-----------------------------------------

    void RK33_compt()
    {
        // compt_t1(cff1.uu_t1, cff1.uu_t0,n_all);
        // uu_t1 = uu_t0 + dt * df;  df_dx(double *g, int i)
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = df_dx_up(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + dt * df[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = df_dx_up(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = u_n0[nb + j] * 0.75 + 0.25 * u_n1[nb + j] + 0.25 * dt * df[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = df_dx_up(f_n2, nb + j);
            u_nn[nb + j] = u_n0[nb + j] * 1.0 / 3.0 + u_n2[nb + j] * 2.0 / 3.0 + 2.0 / 3.0 * dt * df[nb + j];
            j++;
        }
    }

    void SSPRK54_compt_Shu() // Strong stability preserving Runge-Kutta and multistep time discretizations P23
    {
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n0[nb + j] = df_dx_up(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 0.391752226571890 * dt * L_n0[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n1[nb + j] = df_dx_up(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = 0.444370493651235 * u_n0[nb + j] + 0.555629506348765 * u_n1[nb + j] + 0.368410593050371 * dt * L_n1[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n2[nb + j] = df_dx_up(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = 0.620101851488403 * u_n0[nb + j] + 0.379898148511597 * u_n2[nb + j] + 0.251891774271694 * dt * L_n2[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n3[nb + j] = df_dx_up(f_n3, nb + j); //-uu_t1[nb + j];
            u_n4[nb + j] = 0.178079954393132 * u_n0[nb + j] + 0.821920045606868 * u_n3[nb + j] + 0.544974750228521 * dt * L_n3[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }

        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = df_dx_up(f_n4, nb + j);
            u_nn[nb + j] = 0.517231671970585 * u_n2[nb + j] + 0.096059710526147 * u_n3[nb + j] + 0.063692468666290 * dt * L_n3[nb + j] +
                           0.386708617503269 * u_n4[nb + j] + 0.226007483236906 * dt * df[nb + j];

            j++;
        }
    }

    // //-----------------------------------------
    void compt_TDMS_stage_int()
    {
        int i, j;

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = df_dx_up(f_n0, j);
            // L_n4[j] = df_dx1(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxcent(L_n0, j);
            Tu2[j] = Tu1[j], TL_n1[j] = TL_n0[j], TLt_n1[j] = TLt_n0[j];
            Tu1[j] = u_n0[j], TL_n0[j] = L_n0[j], TLt_n0[j] = Lt_n0[j];
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
            ofs << " VARIABLES=x,u,fk " << endl;

            for (int i = nb; i < n1 - nb; i++)
            {

                ofs.precision(10);
                ofs << x[i];
                ofs << "  ";
                ofs << u_nn[i];
                ofs << "  ";
                ofs << f_k[i];
                // for (int j = 0; j < kt; j++)
                // {
                //     ofs << "  ";
                //     ofs.precision(18);
                //     ofs << put_obj[j][i];
                // }

                ofs << endl;
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
    void compt_CFL()
    {
        double *umax;
        umax = max_element(u_n0 + 1, u_n0 + n1 - 1);
        // cout << *umax;
        //  dt = CFL * dx / (*umax);
        // dt=0.02/1.0; /0.5
        dt = CFL * dx;
    }
    ////-----------------------------------------
    int Perimeter()
    {
        cout << " yes ";
        return 2 * t;
    }

    double df_dx_up(double *g, int i)
    {
        double y;
        // 1up
        // y = (*(g + i + 1) - *(g + i)) / dx;

        // 2up
        // y = (*(g + i + 1) - *(g + i - 1)) / dx / 2.0;
        // 5up
        //  y = (*(g + i - 3) / (-30.) + *(g + i - 2) / 4. - *(g + i - 1) + *(g + i) / 3. + *(g + i + 1) / 2. - *(g + i + 2) / 20.) / dx;
        // 7up f'>0
        // y = (*(g + i - 4) * 3. - *(g + i - 3) * 28. + *(g + i - 2) * 126. - *(g + i - 1) * 420. + *(g + i) * 105. + *(g + i + 1) * 252. - *(g + i + 2) * 42. + *(g + i + 3) * 4.) / 420. / dx;
        // 7up f'<0
        y = (-*(g + i - 3) * 4. + *(g + i - 2) * 42. - *(g + i - 1) * 252. - *(g + i) * 105. + *(g + i + 1) * 420. - *(g + i + 2) * 126. + *(g + i + 3) * 28. - *(g + i + 4) * 3.) / 420. / dx;

        // 8up  f'<0
        // y = (-*(g + i - 3) * 5. + *(g + i - 2) * 60. - *(g + i - 1) * 420. - *(g + i) * 378 + *(g + i + 1) * 1050. - *(g + i + 2) * 420. + *(g + i + 3) * 140. - *(g + i + 4) * 30. + *(g + i + 5) * 3.) / 840. / dx;
        // 8up  f'>0
        // y = -(-*(g + i + 3) * 5. + *(g + i + 2) * 60. - *(g + i + 1) * 420. - *(g + i) * 378 + *(g + i - 1) * 1050. - *(g + i - 2) * 420. + *(g + i - 3) * 140. - *(g + i - 4) * 30. + *(g + i - 5) * 3.) / 840. / dx;

        // 9up f'<0
        // y = (*(g + i - 4) * 5. - *(g + i - 3) * 60. + *(g + i - 2) * 360. - *(g + i - 1) * 1680. - *(g + i) * 504. + *(g + i + 1) * 2520. - *(g + i + 2) * 840. + *(g + i + 3) * 240. - *(g + i + 4) * 45. + *(g + i + 5) * 4.) / 2520. / dx;

        // 9up  f'>0
        // y = -(*(g + i + 4) - *(g + i + 3) * 12. + *(g + i + 2) * 72. - *(g + i + 1) * 336. - *(g + i) * 100.8 + *(g + i - 1) * 504. - *(g + i - 2) * 168. + *(g + i - 3) * 48. - *(g + i - 4) * 9. + *(g + i - 5) * 0.8) / 504. / dx;

        // 13 order
        // y = (*(g + i - 6) * 1.803751803752541E-004 - *(g + i - 5) * 2.597402597403613E-003 + *(g + i - 4) * 1.785714285714940E-002 -
        //      *(g + i - 3) * 7.936507936510601E-002 + *(g + i - 2) * 0.267857142857228 - *(g + i - 1) * 0.857142857144015 +
        //      *(g + i) * 2.101706824838899E-012 +
        //      *(g + i + 1) * 0.85714285714181 - *(g + i + 2) * 0.267857142857150 + *(g + i + 3) * 7.936507936510066E-002 -
        //      *(g + i + 4) * 1.785714285714922E-002 + *(g + i + 5) * 2.597402597403614E-003 - *(g + i + 6) * 1.803751803752537E-004) /
        //     dx;

        //  1.803751803752541E-004 -2.597402597403613E-003  1.785714285714940E-002
        //  -7.936507936510601E-002  0.267857142857228      -0.857142857144015
        //   2.101706824838899E-012
        //   0.857142857141818      -0.267857142857150  7.936507936510066E-002
        //    -1.785714285714922E-002  2.597402597403614E-003 -1.803751803752537E-004

        return y;
    }
    double df_dx1(double *g, int i)
    {
        double y;
        // 1up
        // y = (*(g + i + 1) - *(g + i)) / dx;
        // 2up
        // y = (*(g + i + 1) - *(g + i - 1)) / dx / 2.0;
        // 5up
        //  y = (*(g + i - 3) / (-30.) + *(g + i - 2) / 4. - *(g + i - 1) + *(g + i) / 3. + *(g + i + 1) / 2. - *(g + i + 2) / 20.) / dx;
        // 7up
        // y = (*(g + i - 4) * 3. - *(g + i - 3) * 28. + *(g + i - 2) * 126. - *(g + i - 1) * 420. + *(g + i) * 105. + *(g + i + 1) * 252. - *(g + i + 2) * 42. + *(g + i + 3) * 4.) / 420. / dx;
        // 8up  f'<0
        y = (-*(g + i - 3) * 5. + *(g + i - 2) * 60. - *(g + i - 1) * 420. - *(g + i) * 378 + *(g + i + 1) * 1050. - *(g + i + 2) * 420. + *(g + i + 3) * 140. - *(g + i + 4) * 30. + *(g + i + 5) * 3.) / 840. / dx;
        // 8up  f'>0
        // y = -(-*(g + i + 3) * 5. + *(g + i + 2) * 60. - *(g + i + 1) * 420. - *(g + i) * 378 + *(g + i - 1) * 1050. - *(g + i - 2) * 420. + *(g + i - 3) * 140. - *(g + i - 4) * 30. + *(g + i - 5) * 3.) / 840. / dx;
        // 9up f'<0
        // y = (*(g + i - 4) * 5. - *(g + i - 3) * 60. + *(g + i - 2) * 360. - *(g + i - 1) * 1680. - *(g + i) * 504. + *(g + i + 1) * 2520. - *(g + i + 2) * 840. + *(g + i + 3) * 240. - *(g + i + 4) * 45. + *(g + i + 5) * 4.) / 2520. / dx;

        return y;
    }

    double df_dx(double *ff, int ii)
    {
        double y, ss1, ss2;
        int i = ii;
        // 1up
        // y = (*(g + i) - *(g + i - 1)) / dx; //u_t+u_x=0
        // y = (*(g + i + 1) - *(g + i)) / dx; // u_t-u_x=0

        /*   weno7
        //     subroutine hh_weno7P(Ka,Kb,v,hh )           ! Ka=-3,  Kb=3
        //  Use OCFD_constants
        //  implicit none
        double S0, S1, S2, S3, S10, S11, S12, S13, S20, S21, S22, S23, S30, S31, S32, S33, tau,
            ax0, ax1, ax2, ax3, am, q0, q1, q2, q3, vz[2];
        double CC0 = 1. / 35., CC1 = 12. / 35., CC2 = 18. / 35., CC3 = 4. / 35.,
               ax11 = -2. / 6., ax12 = 9. / 6., ax13 = -18. / 6., ax14 = 11. / 6.,
               ax21 = 1. / 6., ax23 = 3. / 6., ax24 = 2. / 6.,
               ax31 = -2. / 6., ax32 = -3. / 6., ax34 = -1. / 6.,
               ax41 = -11. / 6., ax42 = 18. / 6., ax43 = -9. / 6., ax44 = 2. / 6.,
               b12 = 4., b13 = -5., b14 = 2., b22 = -2.,
               b41 = 2., b42 = -5., b43 = 4., c12 = 3.,
               d12 = 13. / 12., d13 = 1043. / 960., d14 = 1. / 12.;
        double e11 = -3. / 12., e12 = 13. / 12., e13 = -23. / 12., e14 = 25. / 12.,
               e21 = 1. / 12., e22 = -5. / 12., e23 = 13. / 12., e24 = 3. / 12.,
               e31 = -1. / 12., e32 = 7. / 12., e33 = 7. / 12., e34 = -1. / 12.,
               e41 = 3. / 12., e42 = 13. / 12., e43 = -5. / 12., e44 = 1. / 12., ep = 1.e-6; //   !! WENO-JS
        // for (int j = 0; j < 2; j++)
        // {
        //     // ! 7th order WENO scheme
        //     // ! 1  阶导数      *(fz + i - 3)
        //     S10 = ax11 * *(fz + i - 3) + ax12 * *(fz + i - 2) + ax13 * *(fz + i - 1) + ax14 * *(fz + i);
        //     S11 = ax21 * *(fz + i - 2) - *(fz + i - 1) + ax23 * *(fz + i) + ax24 * *(fz + i + 1);
        //     S12 = ax31 * *(fz + i - 1) + ax32 * *(fz + i) + *(fz + i + 1) + ax34 * *(fz + i + 2);
        //     S13 = ax41 * *(fz + i) + ax42 * *(fz + i + 1) + ax43 * *(fz + i + 2) + ax44 * *(fz + i + 3);
        //     //  ! 2 阶导数
        //     S20 = -*(fz + i - 3) + b12 * *(fz + i - 2) + b13 * *(fz + i - 1) + b14 * *(fz + i);
        //     S21 = *(fz + i - 1) + b22 * *(fz + i) + *(fz + i + 1);
        //     S22 = *(fz + i) + b22 * *(fz + i + 1) + *(fz + i + 2);
        //     S23 = b41 * *(fz + i) + b42 * *(fz + i + 1) + b43 * *(fz + i + 2) - *(fz + i + 3);
        //     // ! 3 阶导数
        //     S30 = -*(fz + i - 3) + c12 * (*(fz + i - 2) - *(fz + i - 1)) + *(fz + i);
        //     S31 = -*(fz + i - 2) + c12 * (*(fz + i - 1) - *(fz + i)) + *(fz + i + 1);
        //     S32 = -*(fz + i - 1) + c12 * (*(fz + i) - *(fz + i + 1)) + *(fz + i + 2);
        //     S33 = -*(fz + i) + c12 * (*(fz + i + 1) - *(fz + i + 2)) + *(fz + i + 3);

        //     S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
        //     S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
        //     S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
        //     S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;

        //     // !-------WENO J-S----------------------
        //     ax0 = CC0 / ((ep + S0) * (ep + S0));
        //     ax1 = CC1 / ((ep + S1) * (ep + S1));
        //     ax2 = CC2 / ((ep + S2) * (ep + S2));
        //     ax3 = CC3 / ((ep + S3) * (ep + S3));
        //     // !-----------------------------------------------
        //     am = ax0 + ax1 + ax2 + ax3;

        //     // !  4阶差分格式的通量
        //     q0 = e11 * *(fz + i - 3) + e12 * *(fz + i - 2) + e13 * *(fz + i - 1) + e14 * *(fz + i);
        //     q1 = e21 * *(fz + i - 2) + e22 * *(fz + i - 1) + e23 * *(fz + i) + e24 * *(fz + i + 1);
        //     q2 = e31 * *(fz + i - 1) + e32 * *(fz + i) + e33 * *(fz + i + 1) + e34 * *(fz + i + 2);
        //     q3 = e41 * *(fz + i) + e42 * *(fz + i + 1) + e43 * *(fz + i + 2) + e44 * *(fz + i + 3);

        //     // !  由4个4阶差分格式组合成1个7阶差分格式
        //     // !     hj(0)=W0*q0+W1*q1+W2*q2+W3*q3;
        //     vz[j] = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;
        //     i--;
        // }
        // ss1 = (vz[0] - vz[1]) / dx;
        // // -------------------------------
        // i = ii; ////////////////////////

        for (int j = 0; j < 2; j++)
        {
            // !      7th order WENO scheme
            // ! 1  阶导数   *(ff + i - 2)
            S10 = ax11 * *(ff + i + 4) + ax12 * *(ff + i + 3) + ax13 * *(ff + i + 2) + ax14 * *(ff + i + 1);
            S11 = ax21 * *(ff + i + 3) - *(ff + i + 2) + ax23 * *(ff + i + 1) + ax24 * *(ff + i);
            S12 = ax31 * *(ff + i + 2) + ax32 * *(ff + i + 1) + *(ff + i) + ax34 * *(ff + i - 1);
            S13 = ax41 * *(ff + i + 1) + ax42 * *(ff + i) + ax43 * *(ff + i - 1) + ax44 * *(ff + i - 2);
            // ! 2 阶导数
            S20 = -*(ff + i + 4) + b12 * *(ff + i + 3) + b13 * *(ff + i + 2) + b14 * *(ff + i + 1);
            S21 = *(ff + i + 2) + b22 * *(ff + i + 1) + *(ff + i);
            S22 = *(ff + i + 1) + b22 * *(ff + i) + *(ff + i - 1);
            S23 = b41 * *(ff + i + 1) + b42 * *(ff + i) + b43 * *(ff + i - 1) - *(ff + i - 2);
            // ! 3 阶导数
            S30 = -*(ff + i + 4) + c12 * (*(ff + i + 3) - *(ff + i + 2)) + *(ff + i + 1);
            S31 = -*(ff + i + 3) + c12 * (*(ff + i + 2) - *(ff + i + 1)) + *(ff + i);
            S32 = -*(ff + i + 2) + c12 * (*(ff + i + 1) - *(ff + i)) + *(ff + i - 1);
            S33 = -*(ff + i + 1) + c12 * (*(ff + i) - *(ff + i - 1)) + *(ff + i - 2);

            S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
            S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
            S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
            S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;

            // !-------WENO Z----------------------
            ax0 = CC0 * (1.0 + pow((tau / (S0 + ep)), 2));
            ax1 = CC1 * (1.0 + pow((tau / (S1 + ep)), 2));
            ax2 = CC2 * (1.0 + pow((tau / (S2 + ep)), 2));
            ax3 = CC3 * (1.0 + pow((tau / (S3 + ep)), 2));
            // !-------WENO J-S----------------------
            // ax0 = CC0 / ((ep + S0) * (ep + S0));
            // ax1 = CC1 / ((ep + S1) * (ep + S1));
            // ax2 = CC2 / ((ep + S2) * (ep + S2));
            // ax3 = CC3 / ((ep + S3) * (ep + S3));

            // !-----------------------------------------------

            am = ax0 + ax1 + ax2 + ax3;

            // !  4阶差分格式的通量
            q0 = e11 * *(ff + i + 4) + e12 * *(ff + i + 3) + e13 * *(ff + i + 2) + e14 * *(ff + i + 1);
            q1 = e21 * *(ff + i + 3) + e22 * *(ff + i + 2) + e23 * *(ff + i + 1) + e24 * *(ff + i);
            q2 = e31 * *(ff + i + 2) + e32 * *(ff + i + 1) + e33 * *(ff + i) + e34 * *(ff + i - 1);
            q3 = e41 * *(ff + i + 1) + e42 * *(ff + i) + e43 * *(ff + i - 1) + e44 * *(ff + i - 2);

            // !  由4个4阶差分格式组合成1个7阶差分格式
            vz[j] = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;

            i--;
        }
        ss2 = (vz[0] - vz[1]) / dx;

        y = ss2;
        // -------------------------------
        */

        return y;
    }

    double dfdxcent(double *g, int i)
    {
        double y;
        // 1up
        // y = (*(g + i) - *(g + i - 1)) / dx;
        // y = (*(g + i + 1) - *(g + i)) / dx;
        // y = (  *(g + i - 2) *0.5 - *(g + i - 1)*2.0 + *(g + i) *1.5) / dx;
        // 3up
        //  y = (*(g + i - 3) / (-3.) + *(g + i - 2) *1.5 - *(g + i - 1)/3.0 + *(g + i) *11.0/ 6.0) / dx;
        // 5up
        // y = (*(g + i - 3) / (-30.) + *(g + i - 2) / 4. - *(g + i - 1) + *(g + i) / 3. + *(g + i + 1) / 2. - *(g + i + 2) / 20.) / dx;
        // 7up
        // y = (*(g + i - 4) * 3. - *(g + i - 3) * 28. + *(g + i - 2) * 126. - *(g + i - 1) * 420. + *(g + i) * 105. + *(g + i + 1) * 252. - *(g + i + 2) * 42. + *(g + i + 3) * 4.) / 420. / dx;
        // y = (*(g + i - 5) * (-4.) + *(g + i - 4) * 35. - *(g + i - 3) * 140. + *(g + i - 2) * 350. - *(g + i - 1) * 700. + *(g + i) * 3221. + *(g + i + 1) * 140. - *(g + i + 2) * 10.) / 420. / dx;
        // y = (45. * (*(g + i + 1) - *(g + i - 1)) - 21. * (*(g + i + 2) - *(g + i - 2)) + (*(g + i + 3) - *(g + i - 3))) / 60. / dx;
        // y = -(*(g + i + 4) * 3. - *(g + i + 3) * 28. + *(g + i + 2) * 126. - *(g + i + 1) * 420. + *(g + i) * 105. + *(g + i - 1) * 252. - *(g + i - 2) * 42. + *(g + i - 3) * 4.) / 420. / dx;

        // y = (*(g + i + 1) - *(g + i - 1)) / 2.0 / dx;
        // 3 up
        // y = (*(g + i - 1) / (-3.) + *(g + i) / (-2.) + *(g + i + 1) + *(g + i + 2) / (-6.)) / dx;
        // 4 center
        // y = (*(g + i - 2) / (12.) + *(g + i - 1) * (-8. / 12.) + *(g + i + 1) * (8. / 12.) + *(g + i + 2) / (-12.)) / dx;

        // 8up_center
        // y = (*(g + i - 4) * 3. - *(g + i - 3) * 32. + *(g + i - 2) * 168. - *(g + i - 1) * 672. - *(g + i) * 0.0 + *(g + i + 1) * 672. - *(g + i + 2) * 168. + *(g + i + 3) * 32. - *(g + i + 4) * 3.) / 840. / dx;

        //   6 center
        y = (*(g + i - 3) / (-60.) + *(g + i - 2) * (9. / 60.) + *(g + i - 1) * (-45. / 60.) + *(g + i + 1) * (45. / 60.) + *(g + i + 2) * (-9. / 60.) + *(g + i + 3) / (60.)) / dx;

        return y;
    }

    double dfdx(double *g, int i) // U_t-U_x=0
    {
        double y;
        // 1up
        // y = (*(g + i) - *(g + i - 1)) / dx;
        // 5up
        // y = (   *(g + i - 2) * 3. - *(g + i - 1) * 30. - *(g + i) * 20. + *(g + i + 1) * 60. - *(g + i + 2) * 15. + *(g + i + 3) * 2. ) / 60. / dx;

        // 7up
        // y = (-*(g + i - 3) * 4. + *(g + i - 2) * 42. - *(g + i - 1) * 252. - *(g + i) * 105. + *(g + i + 1) * 420. - *(g + i + 2) * 126. + *(g + i + 3) * 28. - *(g + i + 4) * 3.) / 420. / dx;

        // 8up   f'<0
        y = (-*(g + i - 3) * 5. + *(g + i - 2) * 60. - *(g + i - 1) * 420. - *(g + i) * 378 + *(g + i + 1) * 1050. - *(g + i + 2) * 420. + *(g + i + 3) * 140. - *(g + i + 4) * 30. + *(g + i + 5) * 3.) / 840. / dx;

        // 9up
        // y = (*(g + i - 4) - *(g + i - 3) * 12. + *(g + i - 2) * 72. - *(g + i - 1) * 336. - *(g + i) * 100.8 + *(g + i + 1) * 504. - *(g + i + 2) * 168. + *(g + i + 3) * 48. - *(g + i + 4) * 9. + *(g + i + 5) * 0.8) / 504. / dx;

        return y;
    }

    double dfdxx(double *g, int i)
    {
        double y;
        // 2ord
        // y = (*(g + i + 1) - *(g + i) * 2.0 + *(g + i - 1)) / dx / dx;
        // 6ord
        // y = (*(g + i + 3) * 2. - *(g + i + 2) * 27.0 + *(g + i + 1) * 270.0 - *(g + i) * 490.0 + *(g + i - 1) * 270.0 - *(g + i - 2) * 27.0 + *(g + i - 3) * 2.) / dx / dx / 180.0;
        // 8ord
        y = (-*(g + i - 4) * 63. + *(g + i - 3) * 896. - *(g + i - 2) * 7056.0 + *(g + i - 1) * 56448.0 - *(g + i) * 100450.0 + *(g + i + 1) * 56448.0 - *(g + i + 2) * 7056.0 + *(g + i + 3) * 896. - *(g + i + 4) * 63.) / dx / dx / 35280.0;

        return y;
    }

 
    void TVabs(int m)
    {
        double TV = 0.0;
        TVmax = 0.0;
        for (int i = nb; i < n1 - nb; i++) // 间断  int i = nb; i < n1 - nb; i++
        {
            TV = abs(u_nn[i] - 0.5 * sin(PI * (x[i] - t)) - 0.5);
            if (TV > TVmax)
            {
                TVmax = TV;
            }
        }
        TV = TVmax;
        cout << "\n  CFL: " << CFL << ",  TV: " << TV << "   \n";

        string Title = "error.txt";
        ofstream ofs;
        if (m < 2)
        {
            ofs.open(Title, ios::out);
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
        else
        {
            ofs.open(Title, ios::app);
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
    }

    void TVabs1(int m)
    {
        double TV = 0.0, TV1;
        TVmax = 0.0;
        for (int i = nb; i < n1 - nb; i++) // 间断  int i = nb; i < n1 - nb; i++
        {
            TV1 = abs(u_nn[i] - 0.5 * sin(PI * (x[i] - t)) - 0.5);
            TV = TV + TV1;
        }
        TV = TV / (1. * (n1 - 2 * nb - 1));
        cout << "\n  CFL: " << CFL << ",  TV: " << TV << "   \n";

        string Title = "error.txt";
        ofstream ofs;
        if (m < 2)
        {
            ofs.open(Title, ios::out);
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
        else
        {
            ofs.open(Title, ios::app);
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
    }
};
////-----------------------------------------
void comput_tt(cfda &);
int comput_main(int, double, double, int);
int main()
{
    int n, oder_n = 6, begin_cell = 30; // 40
    double t = 2.0, CFL = 0.5;          // 计算时间总长，CFL数  CFL = 0.66 1.45
    // n = begin_cell + 1;
    for (int m = 1; m < oder_n; m++)
    {
        // CFL = CFL + 0.1;
        n = pow(2, m - 1) * begin_cell + 1;
        comput_main(n, t, CFL, m);
    }
    // cout << " \n All comput time: ";
    // printf("%d ms", clock()); //输出运行所费时间，单位毫秒ms
    // system("pause");

    return 2;
}

int comput_main(int n, double t, double CFL, int m)
{
    //
    int nb = 6; // 边界网格结点数（有限差分）
    int n1 = n + 6 * 2;
    double main_df[n + 6 * 2], u_excat[n + 6 * 2];
    double Ltt_n0[n + 6 * 2], x[n + 6 * 2]; // 声明变长数组

    double L_n0[n + 6 * 2], L_n1[n + 6 * 2], L_n2[n + 6 * 2], L_n3[n + 6 * 2], L_n4[n + 6 * 2];      // f=L
    double Lt_n0[n + 6 * 2], Lt_n1[n + 6 * 2], Lt_n2[n + 6 * 2], Lt_n3[n + 6 * 2], Lt_n4[n + 6 * 2]; // L_t
    double u_n0[n + 6 * 2], u_n1[n + 6 * 2], u_n2[n + 6 * 2], u_n3[n + 6 * 2], u_n4[n + 6 * 2];      // 推进步中间时刻
    double u_n5[n + 6 * 2], u_n6[n + 6 * 2], u_n7[n + 6 * 2], u_n8[n + 6 * 2], u_n9[n + 6 * 2];      // 推进步中间时刻
    double f_n0[n + 6 * 2], f_n1[n + 6 * 2], f_n2[n + 6 * 2], f_n3[n + 6 * 2], f_n4[n + 6 * 2];      // 推进步中间时刻
    double f_n5[n + 6 * 2], f_n6[n + 6 * 2], f_n7[n + 6 * 2], f_n8[n + 6 * 2], f_k[n + 6 * 2];       // 推进步中间时刻
    double u_nn[n + 6 * 2], f_nn[n + 6 * 2], Lttx_n0[n + 6 * 2];                                     // 下一时刻（n+1)

    double TL_n0[n + 6 * 2], TL_n1[n + 6 * 2], TL_n2[n + 6 * 2], TL_n3[n + 6 * 2], TL_n4[n + 6 * 2];
    double TLt_n0[n + 6 * 2], TLt_n1[n + 6 * 2], TLt_n2[n + 6 * 2], TLt_n3[n + 6 * 2], TLt_n4[n + 6 * 2];
    double Tu1[n + 6 * 2], Tu2[n + 6 * 2];
    class cfda cff;
    cff.kt = 1;
    cff.A1 = 1 / 3.0;
    cff.A2 = 2 / 3.0;
    cff.t = t;          // 1.0 / cff.PI; //计算总时间
    cff.nx_begin = 0.0; // 网格右端点
    cff.nx_end = 2.0;   // 网格左端点
    cff.CFL = CFL, cff.TVmax = 0.0;

    cff.n = n, cff.nb = nb, cff.n1 = n1;
    cff.Ltt_n0 = Ltt_n0, cff.x = x; // 声明变长数组

    cff.u_n0 = u_n0, cff.u_n1 = u_n1, cff.u_n2 = u_n2, cff.u_n3 = u_n3, cff.u_n4 = u_n4;                     // 推进步中间时刻
    cff.u_n5 = u_n5, cff.u_n6 = u_n6, cff.u_n7 = u_n7, cff.u_n8 = u_n8, cff.u_n9 = u_n9;                     // 推进步中间时刻
    cff.f_n0 = f_n0, cff.f_n1 = f_n1, cff.f_n2 = f_n2, cff.f_n3 = f_n3, cff.f_n4 = f_n4;                     // 推进步中间时刻
    cff.f_n5 = f_n5, cff.f_n6 = f_n6, cff.f_n7 = f_n7, cff.f_n8 = f_n8, cff.f_k = f_k;                       // 推进步中间时刻
                                                                                                             // 声明变长数组
    cff.u_nn = u_nn, cff.f_nn = f_nn, cff.Lttx_n0 = Lttx_n0;                                                 // 下一时刻（n+1)
    cff.L_n0 = L_n0, cff.L_n1 = L_n1, cff.L_n2 = L_n2, cff.L_n3 = L_n3, cff.L_n4 = L_n4;                     // f=L
    cff.Lt_n0 = Lt_n0, cff.Lt_n1 = Lt_n1, cff.Lt_n2 = Lt_n2, cff.Lt_n3 = Lt_n3, cff.Lt_n4 = Lt_n4;           // L_t
                                                                                                             // 存储一时刻（n-1)
    cff.TL_n0 = TL_n0, cff.TL_n1 = TL_n1, cff.TL_n2 = TL_n2, cff.TL_n3 = TL_n3, cff.TL_n4 = TL_n4;           // f=TL
    cff.TLt_n0 = TLt_n0, cff.TLt_n1 = TLt_n1, cff.TLt_n2 = TLt_n2, cff.TLt_n3 = TLt_n3, cff.TLt_n4 = TLt_n4; // TL_t

    cff.Tu1 = Tu1, cff.Tu2 = Tu2;
    cff.df = main_df, cff.u_excat = u_excat,
    cff.intc();
    for (cff.ij = 0; cff.ij < cff.kt; cff.ij++)
    {
        // cff.nt = 1000 * pow(2, cff.ij); //计算时间步数
        comput_tt(cff);
        cff.Store_obj(cff.put_obj[cff.ij], cff.u_nn); // 存储数据
        // cff.Write_obj(cff.ij + 1);                    //存储数据
    }
    cff.Write_obj(-1);
    cff.TVabs1(m);
    return 2;
}

////-----------------------------------------

void comput_tt(cfda &cff1)
{
    int n_all, i3 = 0, tt_flag = 0;
    double th;
    cff1.tt = 0.0;
    n_all = cff1.n1;
    cff1.border();
    cff1.f_eq_u(cff1.f_n0, cff1.u_n0);

    for (int i1 = 0; i1 < 20000000; i1++) // 200 50
    {
        cff1.compt_CFL();

        if (i1 % 2 == 0)
            cff1.dt = cff1.dt * 1.04;
        else
            cff1.dt = cff1.dt * 0.96;

        th = cff1.tt + cff1.dt;
        if (th > cff1.t - 1E-8)
            cff1.dt = cff1.t - cff1.tt, tt_flag = 10;

        if (i1 < 5) // 4, 5
        {
            cff1.compt_TDMS_stage_int();

            // cff1.RK33_compt(); // 3-RK
            // cff1.SSPRK54_compt_Shu(); //// // cff1.RK54LS_compt();

            // cff1.compt23_LLttnn_line();
            cff1.compt34_LLttnn_line();
        }
        else
        {

            // cff1.compt_TDMS_23_ALW_VSS();
            cff1.compt_TDMS_34_ALW_VSS();
        }
        cff1.carr(cff1.u_n0, cff1.u_nn, n_all);
        cff1.f_eq_u(cff1.f_n0, cff1.u_nn);
        // cout << " \n comput time: " << cff1.tt << " comput time n:   " << i1;
        cff1.border();
        // cff1.TVabs_1();
        cff1.tt = cff1.tt + cff1.dt;
        cff1.xt1 = cff1.xt;
        cff1.xt = cff1.dt;
        if (tt_flag > 1)
        {
            cout.precision(18); // 精度为18，一般为6
            cout << " \n comput time: " << cff1.tt << " comput time n:   " << i1 + 1;
            break;
        }
    }
}
