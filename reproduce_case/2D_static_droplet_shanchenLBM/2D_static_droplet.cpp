#include <iostream>
#include <cmath>
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _WIN32
#include <direct.h>
#endif

/* ============================================================
   1. 网格与时间参数
   ============================================================ */
const int nx = 128;           // x方向网格数
const int ny = 128;           // y方向网格数
const int npop = 9;           // D2Q9 离散速度数
const int nsteps = 500;     // 总时间步
const int noutput = 10;     // 输出间隔

/* ============================================================
   2. D2Q9 模型参数
   ============================================================ */
const int cx[] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int cy[] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

const double weights[] = {
    4.0 / 9.0,
    1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
};

double tau = 1.0;     // 弛豫时间
double rhol = 1.95;   // 液相初始密度
double rhog = 0.15;   // 气相初始密度
int radius = 20;      // 初始液滴半径
double g = -5.0;      // 相互作用强度（负值表示吸引）

/* ============================================================
   3. 宏观变量 & 分布函数
   ============================================================ */
double rho[nx * ny];
double u1[nx * ny];
double u2[nx * ny];

double f_mem[nx * ny * npop];
double f2_mem[nx * ny * npop];
double feq[npop];

/* ============================================================
   4. 创建输出文件夹
   ============================================================ */
void createFolder(const std::string& folderName)
{
#ifdef _WIN32
    _mkdir(folderName.c_str());
#else
    mkdir(folderName.c_str(), 0777);
#endif
}

/* ============================================================
   5. 输出函数：保存场变量到文件
   ============================================================ */
void writeFile(const std::string& name, double* f, int n)
{
    std::ofstream fout(name.c_str());
    if (!fout.is_open())
    {
        std::cerr << "Error: cannot open file " << name << std::endl;
        return;
    }

    for (int iY = 0; iY < ny; iY++)
    {
        for (int iX = 0; iX < nx; iX++)
        {
            fout << f[iY * nx + iX] << " ";
        }
        fout << "\n";
    }
    fout << std::endl;
    fout.close();
}

/* ============================================================
   6. 计算气相/液相密度 & 压差（Laplace验证）
   ============================================================ */
void calculateDensities()
{
    double rho_gas = 0.0;
    double rho_liq = 0.0;

    // 左下角区域作为气相
    for (int iX = 0; iX < 10; iX++)
        for (int iY = 0; iY < 10; iY++)
            rho_gas += rho[iY * nx + iX];

    // 中心区域作为液相
    for (int iX = nx / 2 - 5; iX < nx / 2 + 5; iX++)
        for (int iY = ny / 2 - 5; iY < ny / 2 + 5; iY++)
            rho_liq += rho[iY * nx + iX];

    rho_gas /= 100.0;
    rho_liq /= 100.0;

    double p_liq = rho_liq / 3.0 + g / 6.0 * pow((1.0 - exp(-rho_liq)), 2);
    double p_gas = rho_gas / 3.0 + g / 6.0 * pow((1.0 - exp(-rho_gas)), 2);

    std::cout << "Rho_gas = " << rho_gas << "   Rho_liq = " << rho_liq << "\n";
    std::cout << "Delta P = " << p_liq - p_gas << "\n";
}

/* ============================================================
   7. 主函数
   ============================================================ */
int main(int argc, char** argv)
{
    // 命令行输入 tau 和 radius
    if (argc == 3)
    {
        tau = atof(argv[1]);
        radius = atoi(argv[2]);
    }

    double* f = f_mem;
    double* f2 = f2_mem;

    // 按参数创建输出文件夹
    std::stringstream folderStream;
    folderStream << "output_tau" << tau << "_R" << radius;
    std::string outputFolder = folderStream.str();
    createFolder(outputFolder);

    /* ---------------- 初始化 ---------------- */
    for (int i = 0; i < nx * ny; i++)
    {
        int iY = i / nx;
        int iX = i % nx;

        // 圆形液滴初始化（已修正为规范坐标写法）
        if ((iX - nx / 2.0) * (iX - nx / 2.0) +
            (iY - ny / 2.0) * (iY - ny / 2.0) <= radius * radius)
            rho[i] = rhol;
        else
            rho[i] = rhog;

        double dense = rho[i];
        double v1 = 0.0, v2 = 0.0;

        u1[i] = 0.0;
        u2[i] = 0.0;

        double usq = v1 * v1 + v2 * v2;

        feq[0] = 4.0 / 9.0 * dense * (1.0 - 1.5 * usq);
        feq[1] = 1.0 / 9.0 * dense * (1.0 + 3.0 * v1 + 4.5 * v1 * v1 - 1.5 * usq);
        feq[2] = 1.0 / 9.0 * dense * (1.0 + 3.0 * v2 + 4.5 * v2 * v2 - 1.5 * usq);
        feq[3] = 1.0 / 9.0 * dense * (1.0 - 3.0 * v1 + 4.5 * v1 * v1 - 1.5 * usq);
        feq[4] = 1.0 / 9.0 * dense * (1.0 - 3.0 * v2 + 4.5 * v2 * v2 - 1.5 * usq);
        feq[5] = 1.0 / 36.0 * dense * (1.0 + 3.0 * (v1 + v2) + 4.5 * (v1 + v2) * (v1 + v2) - 1.5 * usq);
        feq[6] = 1.0 / 36.0 * dense * (1.0 + 3.0 * (-v1 + v2) + 4.5 * (-v1 + v2) * (-v1 + v2) - 1.5 * usq);
        feq[7] = 1.0 / 36.0 * dense * (1.0 + 3.0 * (-v1 - v2) + 4.5 * (v1 + v2) * (v1 + v2) - 1.5 * usq);
        feq[8] = 1.0 / 36.0 * dense * (1.0 + 3.0 * (v1 - v2) + 4.5 * (v1 - v2) * (v1 - v2) - 1.5 * usq);

        for (int k = 0; k < npop; k++)
        {
            f[9 * i + k] = feq[k];
            f2[9 * i + k] = feq[k];
        }
    }

    time_t start = time(NULL);

    /* ============================================================
       8. 主时间循环
       ============================================================ */
    for (int t = 0; t <= nsteps; t++)
    {
        /* ---------- 计算密度 ---------- */
        for (int i = 0; i < nx * ny; i++)
        {
            rho[i] = 0.0;
            for (int k = 0; k < npop; k++)
                rho[i] += f[9 * i + k];
        }

        /* ---------- 碰撞 + 迁移 ---------- */
        for (int iY = 0; iY < ny; iY++)
        {
            for (int iX = 0; iX < nx; iX++)
            {
                int i = iY * nx + iX;
                double dense = rho[i];

                double fx = 0.0, fy = 0.0;

                /* ----- 计算 Shan-Chen 力 ----- */
                for (int k = 0; k < npop; k++)
                {
                    int iX2 = (iX + cx[k] + nx) % nx;
                    int iY2 = (iY + cy[k] + ny) % ny;

                    double psi_nb = 1.0 - exp(-rho[iY2 * nx + iX2]);

                    fx += weights[k] * cx[k] * psi_nb;
                    fy += weights[k] * cy[k] * psi_nb;
                }

                double psi = 1.0 - exp(-rho[i]);

                fx = -g * psi * fx;
                fy = -g * psi * fy;

                /* ----- 计算速度（含力修正） ----- */
                double v1 = (f[9 * i + 1] - f[9 * i + 3] + f[9 * i + 5] - f[9 * i + 6] - f[9 * i + 7] + f[9 * i + 8]) / dense + fx / (2.0 * dense);
                double v2 = (f[9 * i + 2] - f[9 * i + 4] + f[9 * i + 5] + f[9 * i + 6] - f[9 * i + 7] - f[9 * i + 8]) / dense + fy / (2.0 * dense);

                u1[i] = v1;
                u2[i] = v2;

                /* ----- 外力项（Guo forcing） ----- */
                double fpop[9];
                for (int k = 0; k < npop; k++)
                {
                    fpop[k] = weights[k] * (1.0 - 0.5 / tau) * (
                        (3.0 * (cx[k] - v1) + 9.0 * cx[k] * (cx[k] * v1 + cy[k] * v2)) * fx +
                        (3.0 * (cy[k] - v2) + 9.0 * cy[k] * (cx[k] * v1 + cy[k] * v2)) * fy
                    );
                }

                /* ----- 平衡分布函数 ----- */
                double usq = v1 * v1 + v2 * v2;

                feq[0] = 4.0 / 9.0 * dense * (1.0 - 1.5 * usq);
                feq[1] = 1.0 / 9.0 * dense * (1.0 + 3.0 * v1 + 4.5 * v1 * v1 - 1.5 * usq);
                feq[2] = 1.0 / 9.0 * dense * (1.0 + 3.0 * v2 + 4.5 * v2 * v2 - 1.5 * usq);
                feq[3] = 1.0 / 9.0 * dense * (1.0 - 3.0 * v1 + 4.5 * v1 * v1 - 1.5 * usq);
                feq[4] = 1.0 / 9.0 * dense * (1.0 - 3.0 * v2 + 4.5 * v2 * v2 - 1.5 * usq);
                feq[5] = 1.0 / 36.0 * dense * (1.0 + 3.0 * (v1 + v2) + 4.5 * (v1 + v2) * (v1 + v2) - 1.5 * usq);
                feq[6] = 1.0 / 36.0 * dense * (1.0 + 3.0 * (-v1 + v2) + 4.5 * (-v1 + v2) * (-v1 + v2) - 1.5 * usq);
                feq[7] = 1.0 / 36.0 * dense * (1.0 + 3.0 * (-v1 - v2) + 4.5 * (v1 + v2) * (v1 + v2) - 1.5 * usq);
                feq[8] = 1.0 / 36.0 * dense * (1.0 + 3.0 * (v1 - v2) + 4.5 * (v1 - v2) * (v1 - v2) - 1.5 * usq);

                /* ----- 碰撞 + 迁移 ----- */
                for (int k = 0; k < npop; k++)
                {
                    int iX2 = (iX + cx[k] + nx) % nx;
                    int iY2 = (iY + cy[k] + ny) % ny;

                    f[9 * i + k] += -(f[9 * i + k] - feq[k]) / tau + fpop[k];
                    f2[9 * (iY2 * nx + iX2) + k] = f[9 * i + k];
                }
            }
        }

        std::swap(f, f2);

        /* ---------- 输出 ---------- */
        if (t % noutput == 0)
        {
            std::stringstream name;
            name << t;

            writeFile(outputFolder + "/density" + name.str() + ".dat", rho, nx * ny);
            writeFile(outputFolder + "/velx" + name.str() + ".dat", u1, nx * ny);
            writeFile(outputFolder + "/vely" + name.str() + ".dat", u2, nx * ny);

            std::cout << "Step " << t << "\n";
            
            calculateDensities();
        }
    }

    std::cout << "Total time " << time(NULL) - start << " s\n";
    return 0;
}

//g++ 2D_static_droplet.cpp -o 2D_static_droplet
//./2D_static_droplet 1.0 20