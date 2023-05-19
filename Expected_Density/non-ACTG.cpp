#include <vector>
#include <set>
#include <map>
#include <queue>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>
#include <iostream>
#include <cstdlib>
#include <iterator>
#include <utility>
#include <deque>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <cstdio>
#include <cstring>
#include <stack>
#include <functional>

#include "mod_minimizer.hpp"
#include "window_minimizer.hpp"
#include "syncmer.hpp"


typedef long long ll;
// if you end up using long double, you need to set the floating point notation to fixed, and set the percision to be very high
typedef long double ld;

#define INF 2001001001
#define INF2 2e18
#define MOD 1000000007

#define max3(a, b, c) max(a, max(b, c))
#define min3(a, b, c) min(a, min(b, c))
#define pb push_back 
#define pf push_front
#define f first
#define s second
#define mp make_pair
#define pll pair<ll, ll>
#define pii pair<int, int>
#define tp make_tuple

// first four are north, west, east ,south
int dir1[] = {1, 0, -1, 0, 1, 1, -1, -1};
int dir2[] = {0, 1, 0, -1, 1, -1, 1, -1};

int main() {
	
	std::cout << std::fixed << std::setprecision(8);
	// if you use ld, use the above and don't use string stream

    std::string str;
    
    std::vector<std::string> strs;
    freopen("../non-ACTG.txt", "r", stdin);
    for(int i =0; i < 100; i++){
        std::cin >> str;
        strs.pb(str);
    }

    std::vector<std::vector<double>> mod_min_vec(4, std::vector<double>());
    std::vector<std::vector<double>> wind_min_vec(4, std::vector<double>());
    std::vector<std::vector<double>> sync_vec(4, std::vector<double>());

    uint64_t mods[4] = {109, 128, 1009, 1024};
    unsigned l_winds[4] = {7, 8, 17, 16};

    std::vector<int> kmers(100, 0);

    for(int i =0; i < 100; i++){
        int start = 0;
        while(start + 7 < 1e5){
            bool works = true;
            for(int j =0; j < 8; j++){
                if(strs[i][start+j] == 'N'){
                    works = false;
                    start = start+j;
                    break;
                }
            }
            if(works){
                kmers[i]++;
                
            }
            start++;
        }
        std::cout << kmers[i] << " ";
    }
    std::cout << std::endl;



    for(int i =0; i < 4; i++){
        for(int j =0; j < 100; j++){
            digest::ModMin mm(strs[j], 8, mods[i], 0, 0, 0);
            std::vector<size_t> temp;
            mm.roll_minimizer(100000, temp);
            double am = temp.size();
            am /= kmers[i];
            mod_min_vec[i].pb(am);
        }
    }

    
    for(int i =0; i < 4; i++){
        for(int j =0; j < 100; j++){
            digest::WindowMin wm(strs[j], 8, l_winds[i], 0, 0);
            std::vector<size_t> temp;
            wm.roll_minimizer(100000, temp);
            double am = temp.size();
            am /= kmers[i];

            wind_min_vec[i].pb(am);
        }
    }

    for(int i =0; i < 4; i++){
        for(int j =0; j < 100; j++){
            digest::Syncmer syn(strs[j], 8, l_winds[i], 0, 0);
            std::vector<size_t> temp;
            syn.roll_minimizer(100000, temp);
            double am = temp.size();
            am /= kmers[i];

            sync_vec[i].pb(am);
        }
    }
    freopen("out2.txt", "w", stdout);
    for(int i  = 0; i < 4; i++){
        for(size_t j = 0; j < 100; j++){
            std::cout << mod_min_vec[i][j] << " ";
        }
        std::cout << std::endl;
    }

    for(int i  = 0; i < 4; i++){
        for(size_t j = 0; j < 100; j++){
            std::cout << wind_min_vec[i][j] << " ";
        }
        std::cout << std::endl;
    }

    for(int i  = 0; i < 4; i++){
        for(size_t j = 0; j < 100; j++){
            std::cout << sync_vec[i][j] << " ";
        }
        std::cout << std::endl;
    }

    


	return 0;
}