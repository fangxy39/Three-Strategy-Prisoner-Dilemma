#include <iostream>
#include <fstream>
#include <deque>
#include <cassert>
#include <vector>
#include <cstring>
#include <array>
#include <algorithm>
#include <random>
#include <cmath>
using namespace std;

const unsigned char D = 0; //背叛
const unsigned char C = 1; //合作
const unsigned char B = -1; //漠然
const double threshold = 3; //判定蒙特卡洛结束的阈值

//const double WMIN = 0.1;
//const double WMAX = 1.0;

const double K = 0.1;
//double M = 5;

double R = 1;
double S = 0;
double T = 1.5;
double P = 0;
double sig = 0.3;

inline unsigned char strategy(double per, double nei, double W, double V, double a, double b, double p1, double p2) {
    if (per >= nei){
        if (W == 0){
            if (p1 <= a) return C;
            else return B;
        }
        else{
            if (p1 <= a) return D;
            else return B;
        }
    }
    else {
        if (V == 0){
            if (p2 <= b) return C;
            else return B;
        }
        else{
            if (p2 <= b) return D;
            else return B;
        }
    }
}

inline double payoff(unsigned char a, unsigned char b) {
    if (a == D && b == D) 
        return P;
    if (a == D && b == C) 
        return T;
    if (a == C && b == D) 
        return S;
    if (a == C && b == C)
        return R;
    if (a == B || b == B)
        return sig;
    return 0;
}

inline double E_pay(double per, double nei, double W, double V) {
    if (per >= nei){
        if (W == 0)
            return 10;
        else
            return 11;
    }
    else {
        if (V == 0)
            return 0;
        else
            return 1;
    }
}

inline int round_player(int i, int N) { //边界
    if (i >= N) return  i - N;
    if (i < 0) return i + N;
    return i;
}

int evolution(int L = 200) {
    int N = L * L;
    random_device rd;
    mt19937 g(rd());
    uniform_int_distribution <> dir(0, 3);
    uniform_int_distribution <> uuid(0, 9999);
    uniform_int_distribution <> indi(0, N-1);
    uniform_real_distribution <double> r1(0.0, 1.0);
    uniform_real_distribution <double> r2(0.0, 1.0);
    uniform_real_distribution <double> a(0.0, 1.0);
    uniform_real_distribution <double> b(0.0, 1.0);
    uniform_real_distribution <double> p1(0.0, 1.0);
    uniform_real_distribution <double> p2(0.0, 1.0);
    
    int round_uuid = uuid(g);
    assert(N % 2 == 0);

    vector <double> payoffs,E_payoffs, _strats[4], _nstrats[4];
    vector <unsigned char> checkpoint;
    unsigned char decision[40000][4];
    // vector <deque<double>> mem;
    payoffs.resize(N, 0);
    E_payoffs.resize(N, 0);
    for (int i = 0; i < 4; i++){
        _strats[i].resize(N, 0);
        _nstrats[i].resize(N, 0);
    }
    // mem.resize(N);
    checkpoint.resize(N, 0);

    const int direction[4] = { 1, -1, L, -L };

    // 50% , 50% 
    auto &stratsW = _strats[0];
    for (int i = 0; i < N / 2; i++) stratsW[i] = 1;
    for (int i = N / 2; i < N; i++) stratsW[i] = 0;
    shuffle(stratsW.begin(), stratsW.end(), g);

    auto &stratsa = _strats[1];
    for (int i = 0; i < N; i++) stratsa[i] = a(g);

    auto &stratsV = _strats[2];
    for (int i = 0; i < N / 2; i++) stratsV[i] = 1;
    for (int i = N / 2; i < N; i++) stratsV[i] = 0;
    shuffle(stratsV.begin(), stratsV.end(), g);

    auto &stratsb = _strats[3];
    for (int i = 0; i < N; i++) stratsb[i] = b(g);

    //cout << "var data_" << M << " = { " << endl;
    //cout << "epoch: " << MAX_ITER << "," << endl;
    int P1_cnt = 0;
    int P2_cnt = 0;
    int P3_cnt = 0;
    int P4_cnt = 0;
    vector <int> steady;
    //int steady= threshold + 1;
    steady.resize(20000, 0);
    fill(steady.begin(), steady.end(), 0);
    int count = -1;

    ofstream outfile("output_E.txt");
    if (!outfile.is_open()){
        cerr << "无法打开文件output_E.txt" << endl;
        return 1;
    }
    outfile << "步数\tP1\tP2\tP3\tP4\t合作\t背叛\t孤独\t平均收益\t平均期望收益\t更新人数"<< endl;
    fill(payoffs.begin(), payoffs.end(), 0); //将 payoffs 向量的元素全部置为 0
    fill(E_payoffs.begin(), E_payoffs.end(), 0);
    int check = 0;
    //蒙特卡洛循环
    while(true){
        count++;
        //开始一个蒙特卡洛步
        for (int _iter = 0; _iter < N; _iter++) {
            auto &stratsW = _strats[0];
            auto &stratsa = _strats[1];
            auto &stratsV = _strats[2];
            auto &stratsb = _strats[3];
            auto &next_stratsW = _nstrats[0];
            auto &next_stratsa = _nstrats[1];
            auto &next_stratsV = _nstrats[2];
            auto &next_stratsb = _nstrats[3];

            // 随机选择一个人和所有邻居博弈
            int per = indi(g);
            for (int d = 0; d < 4; d++) {
                const int nei = round_player(direction[d] + per, N);
                unsigned char temp = 0;
                if(d==1 || d == 3){
                    decision[per][d] = strategy(payoffs[per], payoffs[nei], stratsW[per], stratsV[per], stratsa[per], stratsb[per], p1(g), p2(g));
                    decision[nei][d-1] = strategy(payoffs[nei], payoffs[per], stratsW[nei], stratsV[nei], stratsa[nei], stratsb[nei], p1(g), p2(g));
                    temp = decision[nei][d-1];
                }
                if(d==0 || d == 2){
                    decision[per][d] = strategy(payoffs[per], payoffs[nei], stratsW[per], stratsV[per], stratsa[per], stratsb[per], p1(g), p2(g));
                    decision[nei][d+1] = strategy(payoffs[nei], payoffs[per], stratsW[nei], stratsV[nei], stratsa[nei], stratsb[nei], p1(g), p2(g));
                    temp = decision[nei][d+1];
                }
                // 计算收益
                payoffs[per] += payoff(decision[per][d], temp);
                
                //比较强弱用真实而不是期望
                int tem = E_pay(payoffs[per],payoffs[nei],stratsW[per], stratsV[per]);
                
                if (tem == 10){
                    if (stratsV[nei] == 0){
                        double p_R = stratsa[per] * stratsb[nei];
                        E_payoffs[per] += R * p_R + sig *(1-p_R);
                    }
                    else{
                        double p_S = stratsa[per] * stratsb[nei];
                        E_payoffs[per] += S * p_S + sig *(1-p_S);
                    }
                }
                if (tem == 11){
                    if (stratsV[nei] == 0){
                        double p_T = stratsa[per] * stratsb[nei];
                        E_payoffs[per] += T * p_T + sig *(1-p_T);
                    }
                    else{
                        double p_P = stratsa[per] * stratsb[nei];
                        E_payoffs[per] += P * p_P + sig *(1-p_P);
                    }
                }
                if (tem == 0){
                    if (stratsW[nei] == 0){
                        double p_R = stratsb[per] * stratsa[nei];
                        E_payoffs[per] += R * p_R + sig *(1-p_R);
                    }
                    else{
                        double p_S = stratsb[per] * stratsa[nei];
                        E_payoffs[per] += S * p_S + sig *(1-p_S);
                    }
                }
                if (tem == 1){
                    if (stratsW[nei] == 0){
                        double p_T = stratsa[per] * stratsb[nei];
                        E_payoffs[per] += T * p_T + sig *(1-p_T);
                    }
                    else{
                        double p_P = stratsa[per] * stratsb[nei];
                        E_payoffs[per] += P * p_P + sig *(1-p_P);
                    }
                }
            }
            E_payoffs[per] = E_payoffs[per]/4;
            payoffs[per] = payoffs[per]/4;
            //cout << per << " " << payoffs[per] << endl;
            // 情感模仿更新
            int nx = 0;
            /*
            for (auto &&q : mem[per])
                if (q == strats[per]) ++nx;
            double wx = WMAX - (WMAX - WMIN) * nx / M;
            */
            double wx = 1;
            const int d = dir(g);
            const int other_player = round_player(direction[d] + per, N);
            //
            
            double Dp = E_payoffs[per] - E_payoffs[other_player];
            //double Dp = payoffs[per] - payoffs[other_player];
            double q = wx / (1.0 + exp(Dp / K));
            double r_1 = r1(g);
            double r_2 = r2(g);
            if (r_1 <= q && r_2 > q) {
                stratsW[per] = stratsW[other_player];
                stratsa[per] = stratsa[other_player];
            } 
            else if (r_1 > q && r_2 <= q) {
                stratsV[per] = stratsV[other_player];
                stratsb[per] = stratsb[other_player];
            } 
            else if (r_1 <= q && r_2 <= q) {
                stratsW[per] = stratsW[other_player];
                stratsa[per] = stratsa[other_player];
                stratsV[per] = stratsV[other_player];
                stratsb[per] = stratsb[other_player];
            } 
            
            /*
            mem[p].push_back(strats[p]);
            if (mem[p].size() > M) {
                mem[p].pop_front();
            }
            */
            /*
            if (_iter % 10 == 0) {
                cout << "\"epoch" << _iter << "\": [" << endl;
                for (int i = 0; i < N; i++) {
                    if (checkpoint[i] != strats[i]) cout << i << ", ";
                }
                cout << "], " << endl;
                checkpoint = strats;
            }
            */
        }
        // 计算情感类型分布  
        int P1 = P1_cnt;
        int P2 = P2_cnt;
        int P3 = P3_cnt;
        int P4 = P4_cnt;  
        P1_cnt = 0;
        P2_cnt = 0;
        P3_cnt = 0;
        P4_cnt = 0;      
        for (int i = 0; i < N; i ++){
            if (stratsW[i] == 1 && stratsV[i] == 1) ++P1_cnt;
            if (stratsW[i] == 0 && stratsV[i] == 1) ++P2_cnt;
            if (stratsW[i] == 1 && stratsV[i] == 0) ++P3_cnt;
            if (stratsW[i] == 0 && stratsV[i] == 0) ++P4_cnt;
        }    
        steady[count] = abs(P1 - P1_cnt)+abs(P2 - P2_cnt)+abs(P3 - P3_cnt)+abs(P4 - P4_cnt);
        if(count % 100 ==0){
            // 计算策略类型和平均收益
            double pay = 0;
            double E_pay = 0;
            int co = 0;
            int de = 0;
            int lo = 0;
            for (int i = 0; i < N; i ++){
                pay += payoffs[i];
                E_pay += E_payoffs[i];
                for (int d = 0; d < 4; d++) {
                    if (decision[i][d] == C) ++co;
                    if (decision[i][d] == D) ++de;
                    if (decision[i][d] == B) ++lo;
                }
            }    
            pay = pay/N;
            E_pay = E_pay/N;
            // 写入数据
            outfile << count << '\t' << P1_cnt << '\t' << P2_cnt << '\t' << P3_cnt << '\t' << P4_cnt << '\t' << co << '\t' << de << '\t' << lo <<  '\t' << pay <<  '\t' << E_pay << '\t' << steady[count] << endl;
        }
        if (count >= 5){
            check = 0;
            for (int k =0; k <5; k ++){
                check += steady[count-k];
            }
            if (check <= 0){
                cout <<  check << endl;
                break;
            }
        }
    }
    cout << "第" << count << "个蒙特卡洛步结束" << endl;
    // 关闭文件
    outfile.close();
    cout << "数据已写入output_E.txt文件" << endl;
    cout << P1_cnt << " " << P2_cnt << " " << P3_cnt << " " << P4_cnt << endl;
    return 0;
}

int main() {
    return evolution();
}