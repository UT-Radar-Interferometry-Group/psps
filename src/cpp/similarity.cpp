#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

double median(std::vector<double> &v,
              const std::size_t counter){
    std::size_t n = counter/2;
    std::nth_element(v.begin(),v.begin()+n,v.begin()+counter);
    return v[n];
}

std::vector<std::vector<int>> scan_array(
        const unsigned int rdmin,
        const unsigned int rdmax){
    int x,y,p,flag;
    std::vector<std::vector<bool>> visited(
            rdmax,
            std::vector<bool> (rdmax,false));
    std::vector<std::vector<int>> indices;
    visited[0][0] = true;
    for (int r = 1; r < (int)rdmax; ++r){
        x = r;
        y = 0;
        p = 1-r;
        if (r > (int)rdmin){
            indices.push_back(std::vector<int> {r,0});
            indices.push_back(std::vector<int> {-r,0});
            indices.push_back(std::vector<int> {0,r});
            indices.push_back(std::vector<int> {0,-r});
        }
        visited[r][0] = true;
        visited[0][r] = true;
        flag = 0;
        while (x>y){
            // do not need to fill holes
            if (flag == 0){
                y++;
                if (p <= 0){
                    // Mid-point is inside or on the perimeter
                    p += 2*y + 1;
                }else{
                    // Mid-point is outside the perimeter
                    x--;
                    p += 2*y - 2*x + 1;
                }
            }else{
                flag--;
            }

            // All the perimeter points have already been visited
            if (x<y){
                break;
            }
            while (!visited[x-1][y]){
                x--;
                flag++;
            }
            visited[x][y]=true;
            visited[y][x]=true;
            if (r > (int)rdmin){
                indices.push_back(std::vector<int> {x,y});
                indices.push_back(std::vector<int> {-x,-y});
                indices.push_back(std::vector<int> {x,-y});
                indices.push_back(std::vector<int> {-x,y});
                if (x != y){
                    indices.push_back(std::vector<int> {y,x});
                    indices.push_back(std::vector<int> {-y,-x});
                    indices.push_back(std::vector<int> {y,-x});
                    indices.push_back(std::vector<int> {-y,x});
                }
            }
            if (flag > 0){
                x++;
            }
        }
    }
    return indices;
}

std::vector<std::vector<double>> median_similarity(
        const std::vector<std::vector<std::vector<double>>> &stack,
        const std::vector<std::vector<bool>> &ps,
        const unsigned int N,
        const unsigned int rdmin,
        const unsigned int rdmax){
    unsigned int nifg, nrow, ncol, nindices;
    nifg = stack.size();
    nrow = ps.size();
    ncol = ps[0].size();
    std::vector<std::vector<double>> sim_median(
            nrow,std::vector<double>(ncol,0.));
    std::vector<std::vector<int>> indices;
    indices = scan_array(rdmin,rdmax);
    nindices = indices.size();
#pragma omp parallel for shared(stack,ps,sim_median,indices,nifg,nindices)
    for (int r0 = 0; r0 < (int)nrow; ++r0){
        for (int c0 = 0; c0 < (int)ncol; ++c0){
            if (!ps[r0][c0]){
                continue;
            }
            unsigned int counter = 0;
            int r,c; 
            std::vector<double> simvec(N,0.);
            for (std::size_t i = 0; i < nindices; ++i){
                r = r0 + indices[i][0];
                c = c0 + indices[i][1];
                if ((r >= 0) && (r < (int)nrow) && (c >= 0) && (c < (int)ncol) && ps[r][c]){
                    double sim = 0.;
                    for (std::size_t j = 0; j < nifg; ++j){
                        sim += std::cos(stack[j][r0][c0] - stack[j][r][c]);
                    }
                    simvec[counter++] = sim/nifg;
                    if (counter >= N){
                        break;
                    }
                }
            }
            sim_median[r0][c0] = median(simvec,counter);
        }
    }
    return sim_median;
}

std::vector<std::vector<double>> max_similarity(
        const std::vector<std::vector<std::vector<double>>> &stack,
        const std::vector<std::vector<bool>> &ps0,
        const double sim_th,
        const unsigned int N,
        const unsigned int rdmin,
        const unsigned int rdmax,
        const unsigned int maxiter,
        const bool nonps_cal_flag){
    unsigned int nifg, nrow, ncol, nindices;
    nifg = stack.size();
    nrow = ps0.size();
    ncol = ps0[0].size();
    std::vector<std::vector<double>> sim_max(
            nrow,std::vector<double>(ncol,0.));
    std::vector<std::vector<bool>> ps_prev(
            nrow,std::vector<bool>(ncol,false));
    std::vector<std::vector<bool>> ps(ps0);
    std::vector<std::vector<int>> indices;
    indices = scan_array(rdmin,rdmax);
    nindices = indices.size();
    std::cout << nrow << " " << ncol << std::endl;
    std::vector<std::vector<int>> rd2(nrow,
            std::vector<int>(ncol,(int)(rdmax*rdmax)));

    for (std::size_t it = 0; it < maxiter; ++it){
        std::cout << "calculating max similarity iterate " << it << std::endl;
        if (it > 0){
#pragma omp parallel for shared(ps,ps_prev,sim_max)
            for (std::size_t i = 0; i < nrow; ++i){
                for (std::size_t j = 0; j < ncol; ++j){
                    ps_prev[i][j] = ps[i][j];
                    ps[i][j] = (sim_max[i][j] >= sim_th);
                }
            }
        }else{
#pragma omp parallel for shared(stack,ps,sim_max,rd2,indices,nifg,nindices)
            for (int r0 = 0; r0 < (int)nrow; ++r0){
                for (int c0 = 0; c0 < (int)ncol; ++c0){
                    if (ps[r0][c0] & !nonps_cal_flag | !ps[r0][c0] & nonps_cal_flag){
                        continue;
                    }
                    int r,c;
                    unsigned int npts = 0;
                    for (std::size_t i = 0; i < nindices; ++i){
                        r = r0 + indices[i][0];
                        c = c0 + indices[i][1];
                        if ((r >= 0) && (r < (int)nrow) && (c >= 0) && (c < (int)ncol)){
                            if (!ps[r][c] & !nonps_cal_flag){
                                continue;
                            }
                            ++npts;
                            double sim = 0.;
                            for (std::size_t j = 0; j < nifg; ++j){
                                sim += std::cos(stack[j][r0][c0] - stack[j][r][c]);
                            }
                            sim /= nifg;
                            sim_max[r0][c0] = sim>sim_max[r0][c0] ? sim : sim_max[r0][c0];
                        }
                        if (npts >= N){
                            //rd2[r0][c0] = (r-r0)*(r-r0)+(c-c0)*(c-c0);
                            break;
                        }
                    }
                }
            }
        }

        if (nonps_cal_flag){
            break;
        }

        unsigned int newps_counter = 0;
#pragma omp parallel for shared(stack,ps,ps_prev,sim_max,rd2,indices,nifg,nindices,newps_counter)
        for (int r0 = 0; r0 < (int)nrow; ++r0){
            for (int c0 = 0; c0 < (int)ncol; ++c0){
                if (ps_prev[r0][c0] || !ps[r0][c0]){
                    continue;
                }
                ++newps_counter;
                int r,c;
                for (std::size_t i = 0; i < nindices; ++i){
                    r = r0 + indices[i][0];
                    c = c0 + indices[i][1];
                    if ((r >= 0) && (r < (int)nrow) && (c >= 0) && (c < (int)ncol)){
                        //int dist = (r-r0)*(r-r0) + (c-c0)*(c-c0);
                        //if (dist > rd2[r][c]){
                        //    continue;
                        //}
                        double sim = 0.;
                        for (std::size_t j = 0; j < nifg; ++j){
                            sim += std::cos(stack[j][r0][c0] - stack[j][r][c]);
                        }
                        sim /= nifg;
                        sim_max[r][c] = sim>sim_max[r][c] ? sim : sim_max[r][c];
                    }
                }
            }
        } // end loop over all pixels
        // check if there is any new PS in this iteration
        if (newps_counter == 0){
            break;
        }
    } // end loop iteration
    return sim_max;
}

int main(){
    return 0;
}


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
PYBIND11_MODULE(cppsimilarity,m){
    m.doc() = R"pbdoc(
        PSInSAR package
        )pbdoc";

    m.def("max_similarity",&max_similarity,R"pbdoc(
        Calculate maximum phase similarity
        )pbdoc");

    m.def("median_similarity",&median_similarity,R"pbdoc(
        Calculate median phase similarity
        )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
