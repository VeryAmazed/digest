#include <bits/stdc++.h>

#include "digest/mod_minimizer.hpp"
#include "digest/window_minimizer.hpp"
#include "digest/syncmer.hpp"
#include "digest/data_structure.hpp" 

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
    freopen("../tests/density/ACTG.txt", "r", stdin);
    for(int i =0; i < 100; i++){
        std::cin >> str;
        strs.pb(str);
    }

    std::vector<std::vector<double>> mod_min_vec(4, std::vector<double>());
    std::vector<std::vector<double>> wind_min_vec(4, std::vector<double>());
    std::vector<std::vector<double>> sync_vec(4, std::vector<double>());

    uint64_t mods[4] = {109, 128, 1009, 1024};
    unsigned l_winds[4] = {7, 8, 17, 16};

    double kmers = 100000 - 16 + 1;

    for(int i =0; i < 4; i++){
        for(int j =0; j < 100; j++){
            digest::ModMin mm(strs[j], 16, mods[i], 0, 0, digest::MinimizedHashType::CANON);
            std::vector<uint32_t> temp;
            mm.roll_minimizer(100000, temp);
            double am = temp.size();
            am /= kmers;
            mod_min_vec[i].pb(am);
        }
    }

    
    for(int i =0; i < 4; i++){
        for(int j =0; j < 100; j++){
            digest::WindowMin<data_structure::Adaptive> wm(strs[j], 16, l_winds[i], 0, digest::MinimizedHashType::CANON);
            std::vector<uint32_t> temp;
            wm.roll_minimizer(100000, temp);
            double am = temp.size();
            am /= kmers;
            wind_min_vec[i].pb(am);
        }
    }

    for(int i =0; i < 4; i++){
        for(int j =0; j < 100; j++){
            digest::Syncmer<data_structure::Adaptive> syn(strs[j], 16, l_winds[i], 0, digest::MinimizedHashType::CANON);
            std::vector<uint32_t> temp;
            syn.roll_minimizer(100000, temp);
            double am = temp.size();
            am /= kmers;
            sync_vec[i].pb(am);
        }
    }
    freopen("../tests/density/out1.txt", "w", stdout);
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