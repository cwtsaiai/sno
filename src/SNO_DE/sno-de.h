#ifndef __SNO_DE_H_INCLUDED__
#define __SNO_DE_H_INCLUDED__

#include <iostream>
#include <fstream>
#include <limits>
#include <functional>
#include <algorithm> 
#include <random>
#include <filesystem>
#include <iomanip>
#include "test_function.h"

using namespace std;

class sno
{
public:
    using solution = vector<double>;
    using population = vector<solution>;

    struct expected_value
    {
        int I_a;
        int I_b;
        double I_ratio;
        double improvement;
        double region_best_val;
        double final_value;
    };
    using ept_value = expected_value;

    struct region
    {
        vector<int> member_index;
        int best_member;
        ept_value epv_of_region;
    };

    struct net
    {
        population elastic_points;
        vector<double> obj_vals_eps;
        vector<double> improvement_eps;
        vector<region> set_regions;
    };

    sno(int num_runs,
	int num_evals,
	int num_dims,
	int net_size,
	int pop_size_s_init,
	int pop_size_s_end,
	int pop_size_x_init,
	int pop_size_x_end,
	double de_f,
	double de_cr,
	int num_adjustment_max,
	int func_num
	);

    sno();
    virtual ~sno() {}

    void run();

    double evaluate(int dims, solution& candidate_sol);
    void cal_expected_value(net& space_net, const int eval_count);
    void region_search(population& curr_explorer, population& curr_miner, net& space_net, int& eval_count);
    void point_search(population& curr_explorer, population& curr_miner, net& space_net, int& eval_count);
    bool space_net_adjustment(const solution& sol, const double& obj_val,  population& explorer, population& miner, net& space_net, int& eval_count);
    void population_adjustment(const net& space_net, population& curr_explorer, population& curr_miner, int& eval_count);

private:
    void init(population& curr_explorer, population& curr_miner, net& space_net, int& eval_count);
    void region_construct(net& space_net);
    double uniform_rand();
    template <typename T, typename U> vector<int> get_order_index(const vector<T>& value, const U vcomp);
    template <typename T> double euclidean_dis(const T* a, const T* b, const int dims);
    vector< pair<double, int> > roulette_wheel_construct(const vector<int>& ept_rank, const net& space_net, const int eval_count);
    int region_select(const vector< pair<double, int> >& roulette_wheel);
    double boundary_check(const double& value);

private:
    int num_runs;
    int num_evals;
    int num_dims;

    double best_obj_val;
    solution best_sol;
    double avg_obj_val;
    vector<double> avg_obj_val_eval;
    vector<double> run_obj_val;

    int net_size; // size of space net, the size of elastic points will be "net_size * net_size"
    int pop_size_s_init; // initial size of explorer
    int pop_size_s_end; // minimum size of explorer
    int pop_size_x_init; // initial size of miner
    int pop_size_x_end; // minimum size of miner
    vector<double> curr_obj_vals_s; //  objective value of explorer
    vector<double> curr_obj_vals_x; // objective value of miner

    double de_f;
    double de_cr;
    int num_adjustment_max; // for space net adjustment
    // search range
    double bound_min; 
    double bound_max;
    // function information
    int func_num;
    // random number generator
    mt19937_64 rg;
};

// default constructor 
inline sno::sno()
    : num_runs(51),
      num_evals(100000),
      num_dims(10),
      net_size(9),
      pop_size_s_init(190),
      pop_size_s_end(10),
      pop_size_x_init(round(0.1 * 190)),
      pop_size_x_end(round(0.2 * 190)),
      de_f(0.7),
      de_cr(0.5),
      num_adjustment_max(5),
      func_num(1)
{
    cout << "Create the SNO with default constructor." << endl;
    rg.seed(time(0));
}

inline sno::sno(int num_runs,
		int num_evals,
		int num_dims,
		int net_size,
		int pop_size_s_init,
		int pop_size_s_end,
		int pop_size_x_init,
		int pop_size_x_end,
		double de_f,
		double de_cr,
		int num_adjustment_max,
		int func_num
		)
    : num_runs(num_runs),
      num_evals(num_evals),
      num_dims(num_dims),
      net_size(net_size),
      pop_size_s_init(pop_size_s_init),
      pop_size_s_end(pop_size_s_end),
      pop_size_x_init(pop_size_x_init),
      pop_size_x_end(pop_size_x_end),
      de_f(de_f),
      de_cr(de_cr),
      num_adjustment_max(num_adjustment_max),
      func_num(func_num)
{
    cout << "# num_runs: " << num_runs << endl
	 << "# num_evals: " << num_evals << endl
	 << "# num_dims: " << num_dims << endl
	 << "# net_size: " << net_size << endl
	 << "# pop_size_s_init: " << pop_size_s_init << endl
	 << "# pop_size_s_end: " << pop_size_s_end << endl
	 << "# pop_size_x_init: " << pop_size_x_init << endl
	 << "# pop_size_x_end: " << pop_size_x_end << endl
	 << "# F: " << de_f << endl
	 << "# Cr: " << de_cr << endl
	 << "# num_adjustment_max: " << num_adjustment_max << endl
	 << "# func_num: " << func_num << endl
	 << endl;
    rg.seed(time(0));
}

inline double sno::uniform_rand()
{
    uniform_real_distribution<double> u(0.0, 1.0);
    return u(rg);
}

template<typename T, typename U>
vector<int> sno::get_order_index(const vector<T>& value, const U vcomp)
{
    // return the index vector of value after sorting

    vector< pair<T, int> > sort_list;
    for (size_t i = 0; i < value.size(); ++i)
        sort_list.push_back(make_pair(value[i], i));

    sort(sort_list.begin(), sort_list.end(), vcomp);

    vector<int> index;
    for (size_t i = 0; i < sort_list.size(); ++i)
        index.push_back(sort_list[i].second);

    return index;
}

template<typename T>
double sno::euclidean_dis(const T* a, const T* b, const int dims)
{
    double distance = 0.0;
    for (int i = 0; i < dims; ++i)
        distance += (b[i] - a[i]) * (b[i] - a[i]);
    distance = sqrt(distance);
    return distance;
}

inline void sno::run()
{
    cout << scientific << setprecision(10);
    set_search_range(&bound_max, &bound_min, func_num);

    population curr_explorer;
    population curr_miner;
    net space_net;

    avg_obj_val = 0.0;
    avg_obj_val_eval.assign(num_evals, 0.0);
    run_obj_val.assign(num_runs, 0.0);

    for (int r = 0; r < num_runs; r++) {
        int eval_count = 0;
        // initialization.
        init(curr_explorer, curr_miner, space_net, eval_count);
        
        while (eval_count < num_evals) {
            // calculate the expected value of each region.
            cal_expected_value(space_net, eval_count);
            // explorer do the region search.
            region_search(curr_explorer, curr_miner, space_net, eval_count);
            // miner do the point search.
            point_search(curr_explorer, curr_miner, space_net, eval_count);
            // adjust the population size of explorer and miner.
            population_adjustment(space_net, curr_explorer, curr_miner, eval_count);
        }
        avg_obj_val += best_obj_val;
        run_obj_val[r] = best_obj_val;

        cout << "run: " << (r+1) << ", function num: " << func_num
	     << ", dimension: " << num_dims
	     << ", search result: " << best_obj_val << endl;
    }
}

inline void sno::init(population& curr_explorer, population& curr_miner, net& space_net, int& eval_count)
{
    best_sol = solution(num_dims, bound_max);
    best_obj_val = numeric_limits<double>::max();

    // initialization of explorer 
    curr_explorer.assign(pop_size_s_init, solution(num_dims, bound_max));
    curr_obj_vals_s.assign(pop_size_s_init, numeric_limits<double>::max());
    for (size_t i = 0; i < curr_explorer.size(); ++i) {
        for (size_t j = 0; j < curr_explorer[i].size(); ++j) 
            curr_explorer[i][j] = bound_min + (bound_max - bound_min) * uniform_rand();
        
        curr_obj_vals_s[i] = evaluate(num_dims, curr_explorer[i]);
        if (++eval_count > num_evals)
            curr_obj_vals_s[i] = numeric_limits<double>::max();
        
        if (curr_obj_vals_s[i] < best_obj_val) {
            best_obj_val = curr_obj_vals_s[i];
            best_sol = curr_explorer[i];
        }

        if (eval_count <= num_evals)
            avg_obj_val_eval[eval_count - 1] += best_obj_val;
    }
    // initialization of miner 
    curr_miner.assign(pop_size_x_init, solution(num_dims, bound_max));
    curr_obj_vals_x.assign(pop_size_x_init, numeric_limits<double>::max());
    for (size_t i = 0; i < curr_miner.size(); ++i) {
        for (size_t j = 0; j < curr_miner[i].size(); ++j) 
            curr_miner[i][j] = bound_min + (bound_max - bound_min) * uniform_rand();

        curr_obj_vals_x[i] = evaluate(num_dims, curr_miner[i]);
        if (++eval_count > num_evals)
            curr_obj_vals_x[i] = numeric_limits<double>::max();

        if (curr_obj_vals_x[i] < best_obj_val) {
            best_obj_val = curr_obj_vals_x[i];
            best_sol = curr_miner[i];
        }

        if (eval_count <= num_evals)
            avg_obj_val_eval[eval_count - 1] += best_obj_val;
    }
    // initialization of space net 
    int eps_size = net_size * net_size;
    space_net.elastic_points.assign(eps_size, solution(num_dims, bound_max));
    space_net.obj_vals_eps.assign(eps_size, numeric_limits<double>::max());
    space_net.improvement_eps.assign(eps_size, 0.0);
    region_construct(space_net);
    for (size_t i = 0; i < space_net.elastic_points.size(); ++i) {
        for (size_t j = 0; j < space_net.elastic_points[i].size(); ++j)
            space_net.elastic_points[i][j] = bound_min + (bound_max - bound_min) * uniform_rand();

        space_net.obj_vals_eps[i] = evaluate(num_dims, space_net.elastic_points[i]);
        if (++eval_count > num_evals)
            space_net.obj_vals_eps[i] = numeric_limits<double>::max();

        if (space_net.obj_vals_eps[i] < best_obj_val) {
            best_obj_val = space_net.obj_vals_eps[i];
            best_sol = space_net.elastic_points[i];
        }

        if (eval_count <= num_evals)
            avg_obj_val_eval[eval_count - 1] += best_obj_val;
    }    
}

inline double sno::evaluate(int dims, solution& candidate_sol)
{
    double obj_vals = numeric_limits<double>::max();
    obj_vals = cal_test_function(&candidate_sol[0], dims, func_num);

    return obj_vals;
}

inline void sno::region_construct(net& space_net)
{
    int region_size = (net_size - 1) * (net_size - 1);
    // initialize the information of region
    space_net.set_regions.assign(region_size, 
				 region{
				     .member_index = vector<int>(4, -1),
				     .best_member = -1,
				     .epv_of_region = ept_value{
					 .I_a = 1,
					 .I_b = 1,
					 .I_ratio = 1.0,
					 .improvement = 0.0,
					 .region_best_val = numeric_limits<double>::max(),
					 .final_value = 0.0
				     }
				 });
    
    // construct the net according to the index of the elastic points
    int index = 0;
    for (int i = 0; i < region_size; ++i) {
        if ((index + 1) % net_size == 0) 
            ++index;
        space_net.set_regions[i].member_index[0] = index;
        space_net.set_regions[i].member_index[1] = index + 1;
        space_net.set_regions[i].member_index[2] = index + net_size;
        space_net.set_regions[i].member_index[3] = index + net_size + 1;
        ++index;
    }
}

inline void sno::cal_expected_value(net& space_net, const int eval_count)
{
    double I_ratio_max(0.0);
    double I_ratio_min(numeric_limits<double>::max());
    double improvement_max(0.0);
    double improvement_min(numeric_limits<double>::max());
    double rbest_val_max(0.0);
    double rbest_val_min(numeric_limits<double>::max());

    int region_size = (net_size - 1) * (net_size - 1);
    for (int i = 0; i < region_size; ++i) {
        // visited ratio
        space_net.set_regions[i].epv_of_region.I_ratio = 
	    static_cast<double>(space_net.set_regions[i].epv_of_region.I_b) / space_net.set_regions[i].epv_of_region.I_a;
        // region improvement and region best value
        space_net.set_regions[i].epv_of_region.improvement = 0.0;
        space_net.set_regions[i].epv_of_region.region_best_val = numeric_limits<double>::max();
        for (int j = 0; j < 4; ++j) {
            if (space_net.obj_vals_eps[ space_net.set_regions[i].member_index[j] ] < space_net.set_regions[i].epv_of_region.region_best_val) {
                space_net.set_regions[i].epv_of_region.region_best_val = space_net.obj_vals_eps[ space_net.set_regions[i].member_index[j] ];
                space_net.set_regions[i].best_member = space_net.set_regions[i].member_index[j];
            }
            space_net.set_regions[i].epv_of_region.improvement += space_net.improvement_eps[ space_net.set_regions[i].member_index[j] ];
        }

        I_ratio_max = max(space_net.set_regions[i].epv_of_region.I_ratio, I_ratio_max);
        I_ratio_min = min(space_net.set_regions[i].epv_of_region.I_ratio, I_ratio_min);
        improvement_max = max(space_net.set_regions[i].epv_of_region.improvement, improvement_max);
        improvement_min = min(space_net.set_regions[i].epv_of_region.improvement, improvement_min);
        rbest_val_max = max(space_net.set_regions[i].epv_of_region.region_best_val, rbest_val_max);
        rbest_val_min = min(space_net.set_regions[i].epv_of_region.region_best_val, rbest_val_min);
    }

    double w = (1.0 - 2.0) * eval_count / num_evals + 2.0;
    for (int i = 0; i < region_size; ++i) {
        // value will be normalized before getting the sum.
        space_net.set_regions[i].epv_of_region.final_value = 0.0;
        if (I_ratio_max != I_ratio_min) 
            space_net.set_regions[i].epv_of_region.final_value += (space_net.set_regions[i].epv_of_region.I_ratio - I_ratio_min) / (I_ratio_max - I_ratio_min);
        
        if (improvement_max != improvement_min)
            space_net.set_regions[i].epv_of_region.final_value += (space_net.set_regions[i].epv_of_region.improvement - improvement_min) / (improvement_max - improvement_min);
        
        if (rbest_val_max != rbest_val_min)
            space_net.set_regions[i].epv_of_region.final_value += w * (1.0 - (space_net.set_regions[i].epv_of_region.region_best_val - rbest_val_min) / (rbest_val_max - rbest_val_min));
    }
}

inline vector< pair<double, int> > sno::roulette_wheel_construct(const vector<int>& ept_rank, const net& space_net, const int eval_count)
{
    int region_size = (net_size - 1) * (net_size - 1);
    // the size of roulette wheel will change over time
    int curr_be_used = max(static_cast<int>(round(((0.1 - 1.0) * eval_count / num_evals + 1.0) * region_size)), 1);

    double sum = 0.0;
    for (int i = 0; i < curr_be_used; ++i)
        sum += space_net.set_regions[ ept_rank[i] ].epv_of_region.final_value;

    // construct the roulette wheel according to the expected value of region (from high value to low value)
    vector< pair<double, int> > roulette_wheel;
    for (int i = 0; i < curr_be_used; ++i)
        roulette_wheel.push_back(make_pair(space_net.set_regions[ ept_rank[i] ].epv_of_region.final_value / sum, ept_rank[i]));

    return roulette_wheel;
}

inline int sno::region_select(const vector< pair<double, int> >& roulette_wheel)
{
    // select a region from roulette wheel
    double rand = uniform_rand();
    int index = (*(roulette_wheel.end() - 1)).second;
    for (size_t i = 0; i < roulette_wheel.size(); ++i) {
        if (rand < roulette_wheel[i].first) {
            index = roulette_wheel[i].second;
            break;
        }
        else {
            rand -= roulette_wheel[i].first;
        }
    }
    return index;
}

inline double sno::boundary_check(const double& value)
{
    return value > bound_max || value < bound_min ? (bound_max - bound_min) * uniform_rand() + bound_min : value;
}

inline void sno::region_search(population& curr_explorer, population& curr_miner, net& space_net, int& eval_count) 
{
    int eps_size = net_size * net_size;
    int region_size = (net_size - 1) * (net_size - 1);
    for (int i = 0; i < eps_size; ++i)
        space_net.improvement_eps[i] *= 0.7;
    vector<double> ept_value(region_size, 0.0);
    for (int i = 0; i < region_size; ++i)
        ept_value[i] = space_net.set_regions[i].epv_of_region.final_value;
    
    // create a roulette wheel according to expected value
    vector<int> ept_rank = get_order_index(ept_value, [](pair<double, int> x, pair<double, int> y){ return x.first > y.first; });
    vector< pair<double, int> > roulette_wheel = roulette_wheel_construct(ept_rank, space_net, eval_count);
    
    population new_explorer(curr_explorer.size(), solution(num_dims, bound_max));
    vector<double> new_obj_vals_s(curr_explorer.size(), numeric_limits<double>::max());
    for (size_t i = 0; i < curr_explorer.size(); ++i) {
        // select a region
        int region_index = region_select(roulette_wheel);
        ++space_net.set_regions[region_index].epv_of_region.I_a;
        space_net.set_regions[region_index].epv_of_region.I_b = 0.0;

        // select an elastic_point from region 
        int eps_index = -1;
        if (uniform_rand() < ((0.1 - 1.0) * eval_count / num_evals + 1.0)) {
            int eps_r1 = round(uniform_rand() * (4 - 1));
            int eps_r2 = round(uniform_rand() * (4 - 1));
            while (eps_r2 == eps_r1) 
                eps_r2 = round(uniform_rand() * (4 - 1));

            if (space_net.obj_vals_eps[ space_net.set_regions[region_index].member_index[eps_r1] ] > space_net.obj_vals_eps[ space_net.set_regions[region_index].member_index[eps_r2] ])
                eps_index = space_net.set_regions[region_index].member_index[eps_r2];
            else
                eps_index = space_net.set_regions[region_index].member_index[eps_r1];
        }
        else {
            eps_index = space_net.set_regions[region_index].best_member;
        }
       
        // randomly select two index for updating
        size_t ri_1, ri_2;
        do {
            ri_1 = round(uniform_rand() * (curr_explorer.size() - 1));
        } while (ri_1 == i);
        do {
            ri_2 = round(uniform_rand() * (curr_explorer.size() - 1));
        } while (ri_2 == i || ri_2 == ri_1);
        
        // generate a new explorer according to update strategy
        int rand_dim = round(uniform_rand() * (num_dims - 1));
        for (int j = 0; j < num_dims; ++j) {
            if (uniform_rand() < de_cr || j == rand_dim) {
                double r1v = curr_explorer[ri_1][j];
                double r2v = curr_explorer[ri_2][j];
                new_explorer[i][j] = r1v + de_f * (space_net.elastic_points[eps_index][j] - r2v);
                
                new_explorer[i][j] = boundary_check(new_explorer[i][j]);
            }
            else {
                new_explorer[i][j] = curr_explorer[i][j];
            }
        }
        
        // evaluate new explorer
        new_obj_vals_s[i] = evaluate(num_dims, new_explorer[i]);
        if (++eval_count > num_evals)
            new_obj_vals_s[i] = numeric_limits<double>::max();

        if (new_obj_vals_s[i] < best_obj_val) {
            best_obj_val = new_obj_vals_s[i];
            best_sol = new_explorer[i];
        }

        if (eval_count <= num_evals)
            avg_obj_val_eval[eval_count - 1] += best_obj_val;
    }
    
    for (size_t i = 0; i < curr_explorer.size(); ++i) {
        if (new_obj_vals_s[i] < curr_obj_vals_s[i]) {
            curr_explorer[i] = new_explorer[i];
            curr_obj_vals_s[i] = new_obj_vals_s[i];
        }
    }

    // adjust space net according to history of region search.
    population search_history = new_explorer;
    for (size_t i = 0; i < search_history.size(); ++i)
        space_net_adjustment(search_history[i], new_obj_vals_s[i], curr_explorer, curr_miner, space_net, eval_count);

    // update region information.
    for (int i = 0; i < region_size; ++i) {
        ++space_net.set_regions[i].epv_of_region.I_b;
        if (space_net.set_regions[i].epv_of_region.I_b > 1)  
            space_net.set_regions[i].epv_of_region.I_a = 1;
    }
}

inline bool sno::space_net_adjustment(const solution& sol, const double& obj_val, population& explorer, population& miner, net& space_net, int& eval_count)
{
    // sort all elastic points according to distance between $sol$ and each elastic point.
    vector< pair<double, size_t> > distance_order;
    for (size_t i = 0; i < space_net.elastic_points.size(); ++i) {
        double dis = euclidean_dis(&sol[0], &space_net.elastic_points[i][0], num_dims);
        distance_order.push_back(make_pair(dis, i));
    }
    sort(distance_order.begin(), distance_order.end());

    // decide how many elastic points can be adjusted in this adjustment.
    int num_adjustment = max(static_cast<int>(round(static_cast<double>(num_adjustment_max) * eval_count / num_evals)), 1);
    bool return_flag = false;

    for (int i = 0; i < num_adjustment; ++i) {
        // if the objective value of $sol$ worse than its neighbor, the adjustment will be terminated ( because it may imply that $sol$ locate in some worse region )
        if (space_net.obj_vals_eps[ distance_order[i].second ] < obj_val) break;
        // the closest elastic point will be replaced directly, if $sol$ better than it.
        if (i == 0) {
            space_net.improvement_eps[ distance_order[i].second ] += 0.3 * (space_net.obj_vals_eps[ distance_order[i].second ] - obj_val);
            space_net.elastic_points[ distance_order[i].second ] = sol;
            space_net.obj_vals_eps[ distance_order[i].second ] = obj_val;

            return_flag = true;
            continue;
        }
        size_t ri_1, ri_2;
        ri_1 = round(uniform_rand() * (explorer.size() + miner.size() - 1));
        do {
            ri_2 = round(uniform_rand() * (explorer.size() + miner.size() - 1));
        } while (ri_2 == ri_1);

        int rand_dim = round(uniform_rand() * (num_dims - 1));
        solution new_position(num_dims, bound_max);
        for (int j = 0; j < num_dims; ++j) {
            double r1v = (ri_1 >= explorer.size()) ? miner[ri_1 - explorer.size()][j] : explorer[ri_1][j];
            double r2v = (ri_2 >= explorer.size()) ? miner[ri_2 - explorer.size()][j] : explorer[ri_2][j];

            if (uniform_rand() < de_cr || j == rand_dim) {
                new_position[j] = r1v + de_f * (sol[j] - r2v);
                new_position[j] = boundary_check(new_position[j]);
            }
            else {
                new_position[j] = space_net.elastic_points[ distance_order[i].second ][j];
            }
        }

        double new_obj_val = evaluate(num_dims, new_position);
        if (++eval_count > num_evals)
            new_obj_val = numeric_limits<double>::max();

        if (space_net.obj_vals_eps[ distance_order[i].second ] > new_obj_val) {
            space_net.improvement_eps[ distance_order[i].second ] += 0.3 * (space_net.obj_vals_eps[ distance_order[i].second ] - new_obj_val);
            space_net.elastic_points[ distance_order[i].second ] = new_position;
            space_net.obj_vals_eps[ distance_order[i].second ] = new_obj_val;

            if (new_obj_val < best_obj_val) {
                best_obj_val = new_obj_val;
                best_sol = new_position;
            }
        }

        if (eval_count <= num_evals)
            avg_obj_val_eval[eval_count - 1] += best_obj_val;
    }
    
    return return_flag;
}

inline void sno::point_search(population& curr_explorer, population& curr_miner, net& space_net, int& eval_count)
{
    bool sorting_flag = false;
    vector<int> obj_vals_rank;
    for (size_t i = 0; i < curr_miner.size(); ++i) {
        // if sorting_flag == true, it means the space net has been changed, so the rank of eps need to be regenerated
        if (sorting_flag || i == 0) 
            obj_vals_rank = get_order_index(space_net.obj_vals_eps, [](pair<double, int> x, pair<double, int> y) -> bool{ return x.first < y.first; });

        int curr_index = round(uniform_rand() * (curr_miner.size() - 1));
        int ri_1, ri_2;
        do {
            ri_1 = round(uniform_rand() * (curr_miner.size() - 1));
        } while (ri_1 == curr_index);
        do {
            ri_2 = round(uniform_rand() * (curr_miner.size() - 1));
        } while (ri_2 == ri_1 || ri_2 == curr_index);

        solution new_sol(num_dims, bound_max);
        int rand_dim = round(uniform_rand() * (num_dims - 1));
        for (int j = 0; j < num_dims; ++j) {
            if (uniform_rand() < de_cr || j == rand_dim) {
                double r1v = curr_miner[ri_1][j];
                double r2v = curr_miner[ri_2][j];

                new_sol[j] = r1v + de_f * (space_net.elastic_points[obj_vals_rank[0]][j]  - r2v);
                new_sol[j] = boundary_check(new_sol[j]);
            }
            else {
                new_sol[j] = curr_miner[curr_index][j];
            }
        }

        double new_obj_val = evaluate(num_dims, new_sol);
        if (++eval_count > num_evals)
            new_obj_val = numeric_limits<double>::max();

        if (new_obj_val < curr_obj_vals_x[curr_index]) {
            curr_obj_vals_x[curr_index] = new_obj_val;
            curr_miner[curr_index] = new_sol;

            if (new_obj_val < best_obj_val) {
                best_obj_val = new_obj_val;
                best_sol = new_sol;
            }
        }

        if (eval_count <= num_evals)
            avg_obj_val_eval[eval_count - 1] += best_obj_val;

        sorting_flag = space_net_adjustment(new_sol, new_obj_val, curr_explorer, curr_miner, space_net, eval_count);
    }
}

inline void sno::population_adjustment(const net& space_net, population& curr_explorer, population& curr_miner, int& eval_count)
{
    double change_rate = pow(static_cast<double>(eval_count) / num_evals, 0.5);
    // decrease the explorer 
    int new_size_s = (pop_size_s_end - pop_size_s_init) * pow(change_rate, 1.0 - change_rate) + pop_size_s_init;
    while (curr_explorer.size() > static_cast<size_t>(new_size_s)) {
        auto del_index = max_element(curr_obj_vals_s.begin(), curr_obj_vals_s.end()) - curr_obj_vals_s.begin();
        curr_explorer.erase(curr_explorer.begin() + del_index);
        curr_obj_vals_s.erase(curr_obj_vals_s.begin() + del_index);
    }
    
    int new_size_x = (pop_size_x_end - pop_size_x_init) * pow(change_rate, 1.0 - change_rate) + pop_size_x_init;
    vector<int> obj_vals_rank = get_order_index(space_net.obj_vals_eps, [](pair<double, int> x, pair<double, int> y) -> bool{ return x.first < y.first; });
    while (curr_miner.size() < static_cast<size_t>(new_size_x)) {
        solution new_sol(num_dims, bound_max);
        int ri_1, ri_2;
        ri_1 = round(uniform_rand() * (curr_miner.size() - 1));
        do {
            ri_2 = round(uniform_rand() * (curr_miner.size() - 1));
        } while (ri_2 == ri_1);

        // create a new solution with some information of elastic points
        for (size_t i = 0; i < new_sol.size(); ++i) {
            if (uniform_rand() < de_cr) {
                new_sol[i] = curr_miner[ri_1][i] + de_f * (space_net.elastic_points[obj_vals_rank[0]][i] - curr_miner[ri_2][i]);
                new_sol[i] = boundary_check(new_sol[i]);
            }
            else {
                new_sol[i] = space_net.elastic_points[obj_vals_rank[0]][i];
            }
        }

        double new_obj_val = evaluate(num_dims, new_sol);
        if (++eval_count > num_evals)
            new_obj_val = numeric_limits<double>::max();
        
        if (new_obj_val < best_obj_val) {
            best_obj_val = new_obj_val;
            best_sol = new_sol;
        } 

        if (eval_count <= num_evals)
            avg_obj_val_eval[eval_count - 1] += best_obj_val;

        // new solution will be directly put into miner population
        curr_miner.push_back(new_sol);
        curr_obj_vals_x.push_back(new_obj_val);
    }
}

#endif
