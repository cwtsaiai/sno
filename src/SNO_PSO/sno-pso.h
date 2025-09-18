#ifndef __SNO_PSO_H_INCLUDED__
#define __SNO_PSO_H_INCLUDED__

#include <iostream>
#include <fstream>
#include <limits>
#include <functional>
#include <algorithm>
#include <random>
#include <filesystem>
#include <iomanip>
#include <numeric>
#include <cassert>
#include "test_function.h"

using namespace std;

class sno
{
public:
    using solution = vector<double>;
    using population = vector<solution>;

    struct region {
        vector<int> rik;
        int rib;
        int ia;
        int ib;
        double improvement;
        double frib;
        double ei;
    };

    struct net {
        population p;
        vector<double> fp;
        vector<double> improvement_p;
        vector<region> r;
    };

    sno();

    sno(int runs,
	int MaxFES,
	int d,
	int net_size,
	int ns_init,
	int ns_end,
	int nx_init,
	int nx_end,
	double w,
	double c1,
	int na_max,
	int func_num
	);

    virtual ~sno() {}

    void run();

    void Initialization(population& s, population& x, net& space_net, int& eval_count);
    void ExpectedValue(net& space_net, const int eval_count);
    void RegionSearch(population& s, population& x, net& space_net, int& eval_count);
    void PointSearch(population& s, population& x, net& space_net, int& eval_count);
    bool SpaceNetAdjustment(const solution& sol, const double obj_val, population& s, population& x, net& space_net, int& eval_count);
    void PopulationAdjustment(population& s, population& x, const net& space_net, int& eval_count);

private:
    void show_parameters();
    double evaluate(const solution& sol,  int& eval_count) const;
    vector<pair<double, int>> roulette_wheel(const net& space_net, const int eval_count) const;
    int select_region(const vector<pair<double, int>>& rw);
    vector<int> get_ranks(const net& space_net, int num_top_regions, const auto& cmp) const;
    template<typename T> double dist(const vector<T>& a, const vector<T>& b) const;
    double clamp_value(double value);
    double uniform_rand(void);
    int uniform_rand(int last);
    tuple<int, int> uniform_rand_pair(int last, int skip = -1);
    double uniform_rand_in_bound();
    double delta(int eval_count) const;

private:
    const int runs;
    const int MaxFES;
    const int d;
    const int net_size;		// size of space net, the size of elastic points will be "net_size * net_size"
    const int ns_init;		// initial size of explorer
    const int ns_end;		// minimum size of explorer
    const int nx_init;		// initial size of miner
    const int nx_end;		// minimum size of miner
    const double w;
    const double c1;
    const int na_max;		// for space net adjustment
    const int func_num;		// function number

    const int np;		// number of elastic points
    const int nr;		// number of regions

    double fbsf;
    solution bsf;
    double fbsf_avg;
    vector<double> fbsf_avg_FES;
    vector<double> fbsf_run;
    vector<double> fs;		// objective values of explorer
    vector<solution> velocities_s;
    vector<double> fx;		// objective values of miner
    vector<solution> velocities_x;
    double min_bound;		// search range
    double max_bound;
    mt19937_64 rg;		// random number generator
};

// default constructor
inline sno::sno()
    : runs(51),
      MaxFES(100000),
      d(10),
      net_size(9),
      ns_init(190),
      ns_end(10),
      nx_init(round(0.1 * 190)),
      nx_end(round(0.2 * 190)),
      w(0.5),
      c1(1.5),
      na_max(5),
      func_num(1),
      np(net_size * net_size),
      nr((net_size - 1) * (net_size - 1))
{
    show_parameters();
    rg.seed(time(0));
}

inline sno::sno(int runs,
		int MaxFES,
		int d,
		int net_size,
		int ns_init,
		int ns_end,
		int nx_init,
		int nx_end,
		double w,
		double c1,
		int na_max,
		int func_num
		)
    : runs(runs),
      MaxFES(MaxFES),
      d(d),
      net_size(net_size),
      ns_init(ns_init),
      ns_end(ns_end),
      nx_init(nx_init),
      nx_end(nx_end),
      w(w),
      c1(c1),
      na_max(na_max),
      func_num(func_num),
      np(net_size * net_size),
      nr((net_size - 1) * (net_size - 1))
{
    show_parameters();
    rg.seed(time(0));
}

inline void
sno::run()
{
    cout << scientific << setprecision(10);
    cerr << scientific << setprecision(10);

    set_search_range(&max_bound, &min_bound, func_num);

    population s;
    population x;
    net space_net;

    fbsf_avg = 0.0;
    fbsf_avg_FES.assign(MaxFES, 0.0);
    fbsf_run.assign(runs, 0.0);
    for (int r = 0; r < runs; ++r) {
        int eval_count = 0;
        // initialization
        Initialization(s, x, space_net, eval_count);
        while (eval_count < MaxFES) {
            // calculate the expected value of each region
            ExpectedValue(space_net, eval_count);
            // do the region search
            RegionSearch(s, x, space_net, eval_count);
            // do the point search
            PointSearch(s, x, space_net, eval_count);
            // adjust the population size of explorer and miner
            PopulationAdjustment(s, x, space_net, eval_count);
        }
        fbsf_avg += fbsf;
        fbsf_run[r] = fbsf;
        cout << r+1 << " " << fbsf << endl;
    }
    cerr << "# average: " << fbsf_avg / runs << endl;
}

inline void
sno::Initialization(population& s, population& x, net& space_net, int& eval_count)
{
    bsf = solution(d, max_bound);
    fbsf = numeric_limits<double>::max();

    // initialize the explorers
    s.assign(ns_init, solution(d, max_bound));
    fs.assign(ns_init, numeric_limits<double>::max());
    velocities_s.assign(ns_init, solution(d, max_bound));
    for (size_t i = 0; i < s.size(); ++i) {
        for (size_t j = 0; j < s[i].size(); ++j) {
            s[i][j] = uniform_rand_in_bound();
            velocities_s[i][j] = uniform_rand_in_bound();
	}

        fs[i] = evaluate(s[i], eval_count);
        if (eval_count <= MaxFES) {
	    if (fs[i] < fbsf) {
		fbsf = fs[i];
		bsf = s[i];
	    }
            fbsf_avg_FES[eval_count - 1] += fbsf;
	}
    }

    // initialize the miners
    x.assign(nx_init, solution(d, max_bound));
    fx.assign(nx_init, numeric_limits<double>::max());
    velocities_x.assign(nx_init, solution(d, max_bound));
    for (size_t i = 0; i < x.size(); ++i) {
        for (size_t j = 0; j < x[i].size(); ++j) {
            x[i][j] = uniform_rand_in_bound();
            velocities_x[i][j] = uniform_rand_in_bound();
        }

        fx[i] = evaluate(x[i], eval_count);
        if (eval_count <= MaxFES) {
	    if (fx[i] < fbsf) {
		fbsf = fx[i];
		bsf = x[i];
	    }
            fbsf_avg_FES[eval_count - 1] += fbsf;
	}
    }

    // initialize the space net
    space_net.p.assign(np, solution(d, max_bound));
    space_net.fp.assign(np, numeric_limits<double>::max());
    space_net.improvement_p.assign(np, 0.0);
    space_net.r.assign(nr, region {
	    .rik = vector<int>(4, -1),
	    .rib = -1,
	    .ia = 1,
	    .ib = 1,
	    .improvement = 0.0,
	    .frib = numeric_limits<double>::max(),
	    .ei = 0.0
	});

    int n = 0;
    for (int i = 0; i < net_size-1; ++i)
	for (int j = i*net_size; j < (i+1)*net_size-1; ++j)
	    space_net.r[n++].rik = { j, j + 1, j + net_size, j + 1 + net_size };

    for (size_t i = 0; i < space_net.p.size(); ++i) {
        for (size_t j = 0; j < space_net.p[i].size(); ++j)
            space_net.p[i][j] = uniform_rand_in_bound();

        space_net.fp[i] = evaluate(space_net.p[i], eval_count);
        if (eval_count <= MaxFES) {
	    if (space_net.fp[i] < fbsf) {
		fbsf = space_net.fp[i];
		bsf = space_net.p[i];
	    }
            fbsf_avg_FES[eval_count - 1] += fbsf;
	}
    }
}

inline void
sno::ExpectedValue(net& space_net, const int eval_count)
{
    double max_ir = 0.0;
    double min_ir = numeric_limits<double>::max();
    double max_improvement = 0.0;
    double min_improvement = numeric_limits<double>::max();
    double max_frib = 0.0;
    double min_frib = numeric_limits<double>::max();

    for (int i = 0; i < nr; ++i) {
        // visited ratio
        const double ir = static_cast<double>(space_net.r[i].ib) / space_net.r[i].ia;
        // region improvement and region best objective value
        space_net.r[i].improvement = 0.0;
        space_net.r[i].frib = numeric_limits<double>::max();
        for (int j = 0; j < 4; ++j) {
            space_net.r[i].improvement += space_net.improvement_p[space_net.r[i].rik[j]];
            if (space_net.fp[space_net.r[i].rik[j]] < space_net.r[i].frib) {
                space_net.r[i].frib = space_net.fp[space_net.r[i].rik[j]];
                space_net.r[i].rib = space_net.r[i].rik[j];
            }
        }
        max_ir = max(ir, max_ir);
        min_ir = min(ir, min_ir);
        max_improvement = max(space_net.r[i].improvement, max_improvement);
        min_improvement = min(space_net.r[i].improvement, min_improvement);
        max_frib = max(space_net.r[i].frib, max_frib);
        min_frib = min(space_net.r[i].frib, min_frib);
    }

    const double w = 2.0 + delta(eval_count) * (1.0 - 2.0);
    for (int i = 0; i < nr; ++i) {
        space_net.r[i].ei = 0.0;
        const double ir = static_cast<double>(space_net.r[i].ib) / space_net.r[i].ia;
        // all the values are normalized before adding
        if (max_ir != min_ir)
            space_net.r[i].ei += (ir - min_ir) / (max_ir - min_ir);
        if (max_improvement != min_improvement)
            space_net.r[i].ei += (space_net.r[i].improvement - min_improvement) / (max_improvement - min_improvement);
        if (max_frib != min_frib)
            space_net.r[i].ei += w * (1.0 - (space_net.r[i].frib - min_frib) / (max_frib - min_frib));
    }
}

inline void
sno::RegionSearch(population& s, population& x, net& space_net, int& eval_count)
{
    for (int i = 0; i < np; ++i)
        space_net.improvement_p[i] *= 0.7;

    // create a roulette wheel according to the expected values
    const vector<pair<double, int>> rw = roulette_wheel(space_net, eval_count);

    for (size_t i = 0; i < s.size(); ++i) {
        // select a region
        const int n = select_region(rw);
        space_net.r[n].ia += 1;
        space_net.r[n].ib = 0;

        // select an elastic point from the region
        int p_sel = space_net.r[n].rib;
        if (uniform_rand() < 1.0 + delta(eval_count) * (0.1 - 1.0)) {
	    const auto [r1, r2] = uniform_rand_pair(4);
            if (space_net.fp[space_net.r[n].rik[r1]] < space_net.fp[space_net.r[n].rik[r2]])
                p_sel = space_net.r[n].rik[r1];
            else
                p_sel = space_net.r[n].rik[r2];
        }

        for (int j = 0; j < d; ++j) {
            velocities_s[i][j] = w * velocities_s[i][j] + c1 * uniform_rand() * (space_net.p[p_sel][j] - s[i][j]);
            s[i][j] = clamp_value(s[i][j] + velocities_s[i][j]);
        }

        fs[i] = evaluate(s[i], eval_count);
        if (eval_count <= MaxFES) {
	    if (fs[i] < fbsf) {
		fbsf = fs[i];
		bsf = s[i];
	    }
            fbsf_avg_FES[eval_count - 1] += fbsf;
	}
    }

    // adjust the space net according to the history of region search
    population& search_history = s;
    for (size_t i = 0; i < search_history.size(); ++i)
        SpaceNetAdjustment(search_history[i], fs[i], s, x, space_net, eval_count);

    // update region information
    for (int i = 0; i < nr; ++i) {
        space_net.r[i].ib += 1;
        if (space_net.r[i].ib > 1)
            space_net.r[i].ia = 1;
    }
}

inline void
sno::PointSearch(population& s, population& x, net& space_net, int& eval_count)
{
    int top = distance(space_net.fp.begin(), min_element(space_net.fp.begin(), space_net.fp.end()));
    for (size_t i = 0; i < x.size(); ++i) {
	const int phi_x = uniform_rand(x.size());
        for (int j = 0; j < d; ++j) {
            velocities_x[phi_x][j] = w * velocities_x[phi_x][j] + c1 * uniform_rand() * (space_net.p[top][j] - x[phi_x][j]);
            x[phi_x][j] += velocities_x[phi_x][j];
            x[phi_x][j] = clamp_value(x[phi_x][j]);
        }

        fx[phi_x] = evaluate(x[phi_x], eval_count);
        if (eval_count <= MaxFES) {
	    if (fx[phi_x] < fbsf) {
		fbsf = fx[phi_x];
		bsf = x[phi_x];
	    }
            fbsf_avg_FES[eval_count - 1] += fbsf;
	}

        if (SpaceNetAdjustment(x[phi_x], fx[phi_x], s, x, space_net, eval_count))
	    // regenerate the rank of eps if space net adjusted
	    top = distance(space_net.fp.begin(), min_element(space_net.fp.begin(), space_net.fp.end()));
    }
}

inline bool
sno::SpaceNetAdjustment(const solution& sol, const double obj_val, population& s, population& x, net& space_net, int& eval_count)
{
    // partial sort all the elastic points by the distance between sol and elastic point.
    static vector<pair<double, size_t>> distances;
    distances.clear();
    distances.reserve(space_net.p.size());
    for (size_t i = 0; i < space_net.p.size(); ++i)
        distances.push_back(make_pair(dist(sol, space_net.p[i]), i));
    partial_sort(distances.begin(), distances.begin() + na_max, distances.end());

    // decide the number of elastic points to be adjusted.
    const int num_adjusts = max(static_cast<int>(ceil(delta(eval_count) * na_max)), 1);

    bool closest_ep_replaced = false;
    for (int i = 0; i < num_adjusts; ++i) {
	auto& k = distances[i].second;
        // stop the adjustment if sol is worse than its neighbor (because this may imply sol is locates in a bad region.)
        if (space_net.fp[k] < obj_val)
	    break;
        // use sol to replace the closest elastic point if sol is better.
        if (i == 0) {
            space_net.improvement_p[k] += 0.3 * (space_net.fp[k] - obj_val);
            space_net.p[k] = sol;
            space_net.fp[k] = obj_val;
            closest_ep_replaced = true;
            continue;
        }

        solution new_sol(d, max_bound);
        for (int j = 0; j < d; ++j)
            new_sol[j] = clamp_value(space_net.p[k][j] + c1 * uniform_rand() * (sol[j] - space_net.p[k][j]));

        double new_obj_val = evaluate(new_sol, eval_count);
        if (eval_count <= MaxFES) {
            if (new_obj_val < fbsf) {
                fbsf = new_obj_val;
                bsf = new_sol;
            }
            fbsf_avg_FES[eval_count - 1] += fbsf;

	    if (new_obj_val < space_net.fp[k]) {
		space_net.improvement_p[k] += 0.3 * (space_net.fp[k] - new_obj_val);
		space_net.p[k] = new_sol;
		space_net.fp[k] = new_obj_val;
	    }
	}
    }
    return closest_ep_replaced;
}

inline void
sno::PopulationAdjustment(population& s, population& x, const net& space_net, int& eval_count)
{
    const double cr = delta(eval_count);
    
    // decrease the explorers
    // const size_t new_size_s = ns_init + pow(cr, 1.0 - sqrt(cr)) * (ns_end - ns_init);
    const size_t new_size_s = ns_init + pow(sqrt(cr), 1.0 - sqrt(cr)) * (ns_end - ns_init);
    while (s.size() > new_size_s) {
        auto del_index = max_element(fs.begin(), fs.end()) - fs.begin();
        s.erase(s.begin() + del_index);
        fs.erase(fs.begin() + del_index);
        velocities_s.erase(velocities_s.begin() + del_index);
    }

    // increase the miners
    // const size_t new_size_x = nx_init + pow(cr, 1.0 - sqrt(cr)) * (nx_end - nx_init);
    const size_t new_size_x = nx_init + pow(sqrt(cr), 1.0 - sqrt(cr)) * (nx_end - nx_init);
    const int top = distance(space_net.fp.begin(), min_element(space_net.fp.begin(), space_net.fp.end()));
    while (x.size() < new_size_x) {
        // create a new solution with some information taking from elastic points
        solution new_sol(d, max_bound);
        solution new_velocity(d, max_bound);
        const int phi_x = uniform_rand(x.size());
        for (size_t i = 0; i < new_sol.size(); ++i) {
            new_velocity[i] = w * velocities_x[phi_x][i] + c1 * uniform_rand() * (space_net.p[top][i] - x[phi_x][i]);
            new_sol[i] = clamp_value(uniform_rand_in_bound() + new_velocity[i]);
        }

        double new_obj_val = evaluate(new_sol, eval_count);
        if (eval_count <= MaxFES) {
	    if (new_obj_val < fbsf) {
		fbsf = new_obj_val;
		bsf = new_sol;
	    }
            fbsf_avg_FES[eval_count - 1] += fbsf;
	}

        // put the new solution into the miner population
        x.push_back(new_sol);
        fx.push_back(new_obj_val);
        velocities_x.push_back(new_velocity);
    }
}



// helping functions
inline void
sno::show_parameters()
{
    cerr << "# runs: " << runs << endl
	 << "# MaxFES: " << MaxFES << endl
	 << "# d: " << d << endl
	 << "# net_size: " << net_size << endl
	 << "# ns_init: " << ns_init << endl
	 << "# ns_end: " << ns_end << endl
	 << "# nx_init: " << nx_init << endl
	 << "# nx_end: " << nx_end << endl
	 << "# w: " << w << endl
	 << "# cl: " << c1 << endl
	 << "# na_max: " << na_max << endl
	 << "# func_num: " << func_num << endl;
}

inline double
sno::evaluate(const solution& sol, int& eval_count) const
{
    return ++eval_count, cal_test_function(&sol[0], sol.size(), func_num);
}

inline vector<pair<double, int>>
sno::roulette_wheel(const net& space_net, const int eval_count) const
{
    // adapt the size of the roulette wheel
    const int num_top_regions = max(static_cast<int>(round((1.0 + delta(eval_count) * (0.1 - 1.0)) * nr)), 1);

    // get the ranks of top regions
    const vector<int> ranks = get_ranks(space_net, num_top_regions, [](pair<double, int> x, pair<double, int> y) -> bool { return x.first > y.first; });

    // construct the roulette wheel using the expected values of the num_top_regions regions (from high to low)
    double sum = 0.0;
    for (int i = 0; i < num_top_regions; ++i)
        sum += space_net.r[ranks[i]].ei;
    vector<pair<double, int>> rw;
    rw.reserve(num_top_regions);
    for (int i = 0; i < num_top_regions; ++i)
        rw.push_back(make_pair(space_net.r[ranks[i]].ei / sum, ranks[i]));

    return rw;
}

inline int
sno::select_region(const vector<pair<double, int>>& rw)
{
    // select a region from the roulette wheel
    double r = uniform_rand();
    for (size_t i = 0; i < rw.size(); ++i) {
        if (r < rw[i].first)
            return rw[i].second;
        else
            r -= rw[i].first;
    }
    return rw.rbegin()->second;
}

inline vector<int>
sno::get_ranks(const net& space_net, int num_top_regions, const auto& cmp) const
{
    vector<pair<double, int>> vi;
    vi.reserve(space_net.r.size());
    for (size_t i = 0; i < space_net.r.size(); ++i)
        vi.push_back(make_pair(space_net.r[i].ei, i));
    partial_sort(vi.begin(), vi.begin() + num_top_regions, vi.end(), cmp);
    vector<int> ranks;
    ranks.reserve(num_top_regions);
    for (int i = 0; i < num_top_regions; ++i)
        ranks.push_back(vi[i].second);
    return ranks;
}

template<typename T>
inline double
sno::dist(const vector<T>& a, const vector<T>& b) const
{
    assert(a.size() == b.size());
    double v = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        v += (b[i] - a[i]) * (b[i] - a[i]);
    return sqrt(v);
}

inline double
sno::clamp_value(double value)
{
    return value >= min_bound && value <= max_bound ? value : uniform_rand_in_bound();
}

inline double
sno::uniform_rand(void)
{
    return uniform_real_distribution<double>(0.0, 1.0)(rg);
}

inline int
sno::uniform_rand(int last)
{
    return round(uniform_rand() * (last - 1));
}

inline tuple<int, int>
sno::uniform_rand_pair(int last, int skip)
{
    int r1;
    do {
	r1 = uniform_rand(last);
    } while (r1 == skip);
    int r2;
    do {
	r2 = uniform_rand(last);
    } while (r2 == skip || r2 == r1);
    return {r1, r2};
}

inline double
sno::uniform_rand_in_bound()
{
    return min_bound + uniform_rand() * (max_bound - min_bound);
}

inline double
sno::delta(int eval_count) const
{
    return static_cast<double>(eval_count) / MaxFES;
}

#endif
