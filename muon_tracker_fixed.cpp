// muon_tracker_fixed.cpp
// Compile: g++ -O2 -std=c++17 -o muon_tracker_fixed muon_tracker_fixed.cpp
// Run: ./muon_tracker_fixed
//
// Reads hits_0_with_radius.csv and writes tracked_output_top_best_two_iterations_cpp.csv
// Ensures only one best track per top-layer hit (hA_top) across both iterations.

#include <bits/stdc++.h>
using namespace std;

struct Hit {
    int TDCID;
    int CHNLID;
    int eventid;
    int triggerledge;
    double drift_time;
    double adc_time;
    double corr_time;
    double drift_radius; // mm
};

struct SavedHit {
    int track_id;
    int TDCID;
    int CHNLID;
    int eventid;
    int triggerledge;
    double adc_time;
    double drift_time;
    double corr_time;
    double Dt;
    double x;
    double y;
    double drift_radius;
    double residual;
    double a,b,c;
    double chi2ndf;
};

// ---------- geometry ----------
using Point = pair<double,double>;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const string INPUT_CSV = "hit_rad_100k.csv";
    const string OUTPUT_CSV = "tracked_100k.csv";
    const int WINDOW_SIZE = 2000;
    const double CHI2NDF_CUT = 50.0;

    cout << "Loading CSV: " << INPUT_CSV << "\n";
    ifstream fin(INPUT_CSV);
    if (!fin.is_open()) {
        cerr << "Failed to open " << INPUT_CSV << "\n";
        return 1;
    }

    // read header, find indices
    string header;
    if (!getline(fin, header)) {
        cerr << "Empty CSV\n";
        return 1;
    }
    vector<string> cols;
    {
        string cur;
        stringstream ss(header);
        while (getline(ss, cur, ',')) {
            while (!cur.empty() && isspace((unsigned char)cur.back())) cur.pop_back();
            while (!cur.empty() && isspace((unsigned char)cur.front())) cur.erase(cur.begin());
            cols.push_back(cur);
        }
    }
    auto find_col = [&](const string &name)->int {
        for (size_t i=0;i<cols.size();++i){
            string low = cols[i];
            for (auto &c: low) c = tolower((unsigned char)c);
            string n = name;
            for (auto &c: n) c = tolower((unsigned char)c);
            if (low == n) return (int)i;
        }
        for (size_t i=0;i<cols.size();++i){
            string low = cols[i];
            for (auto &c: low) c = tolower((unsigned char)c);
            string n = name;
            for (auto &c: n) c = tolower((unsigned char)c);
            if (low.find(n) != string::npos) return (int)i;
        }
        return -1;
    };

    int idx_TDCID = find_col("TDCID");
    int idx_CHNLID = find_col("CHNLID");
    int idx_eventid = find_col("eventid");
    int idx_triggerledge = find_col("triggerledge");
    int idx_drift = find_col("drift_radius");
    double idx_drift_time = find_col("drift_time");
    double idx_corr_time = find_col("corr_time");
    double idx_adc_time = find_col("adc_time");
    if (idx_drift == -1) idx_drift = find_col("driftradius");
    if (idx_drift == -1) idx_drift = find_col("drift_radius_mm");

    if (idx_TDCID<0 || idx_CHNLID<0 || idx_eventid<0 || idx_triggerledge<0 || idx_drift<0) {
        cerr << "Required columns not found. Found header columns:\n";
        for (auto &c: cols) cerr << c << " | ";
        cerr << "\nNeed TDCID, CHNLID, eventid, triggerledge, drift_radius (or similar)\n";
        return 1;
    }

    // parse CSV into vector<Hit>
    vector<Hit> all_hits;
    string line;
    int line_no = 1;
    while (getline(fin, line)) {
        ++line_no;
        if (line.empty()) continue;
        vector<string> tokens;
        string cur;
        stringstream ss(line);
        while (getline(ss, cur, ',')) tokens.push_back(cur);
        if ((int)tokens.size() < (int)cols.size()) tokens.resize(cols.size());
        try {
            Hit h;
            h.TDCID = stoi(tokens[idx_TDCID]);
            h.CHNLID = stoi(tokens[idx_CHNLID]);
            h.eventid = stoi(tokens[idx_eventid]);
            h.triggerledge = stoi(tokens[idx_triggerledge]);
            h.drift_radius = stod(tokens[idx_drift]);
	    h.drift_time = stod(tokens[idx_drift_time]);
	    h.corr_time = stod(tokens[idx_corr_time]);
	    h.adc_time = stod(tokens[idx_adc_time]);
            all_hits.push_back(h);
        } catch (...) {
            cerr << "Warning: cannot parse line " << line_no << " -> skipping\n";
            continue;
        }
    }
    fin.close();
    cout << "Loaded " << all_hits.size() << " hits\n";

    // geometry setup (cm -> mm)
    map<int, Point> geometry_layer = {
        {0,{1.5,1.5}},{1,{4.5,1.5}},{2,{7.5,1.5}},{3,{10.5,1.5}},
        {4,{13.5,1.5}},{5,{16.5,1.5}},{6,{19.5,1.5}},{7,{22.5,1.5}},
        {8,{3.0,4.1}},{9,{6.0,4.1}},{10,{9.0,4.1}},{11,{12.0,4.1}},
        {12,{15.0,4.1}},{13,{18.0,4.1}},{14,{21.0,4.1}},{15,{24.0,4.1}},
        {16,{1.5,6.7}},{17,{4.5,6.7}},{18,{7.5,6.7}},{19,{10.5,6.7}},
        {20,{13.5,6.7}},{21,{16.5,6.7}},{22,{19.5,6.7}},{23,{22.5,6.7}}
    };
    vector<pair<double,double>> tdc_offsets = {
        {-96.0,0}, {-96.0,34.7}, {-72.0,0}, {-72.0,34.7},
        {-48.0,0}, {-48.0,34.7}, {-24.0,0}, {-24.0,34.7},
        {0.0,0.0}, {0.0,34.7}, {24.0,0.0}, {24.0,34.7},
        {48.0,0.0}, {48.0,34.7}, {72.0,0.0}, {72.0,34.7},
        {96.0,0.0}, {96.0,34.7}
    };
    vector<array<Point,24>> geometry_tdcs;
    geometry_tdcs.resize(tdc_offsets.size());
    for (size_t tid=0; tid<tdc_offsets.size(); ++tid) {
        double dx = tdc_offsets[tid].first;
        double dy = tdc_offsets[tid].second;
        for (int ch=0; ch<24; ++ch) {
            auto p = geometry_layer[ch];
            double xmm = (p.first + dx) * 10.0;
            double ymm = (p.second + dy) * 10.0;
            geometry_tdcs[tid][ch] = {xmm,ymm};
        }
    }

    vector<pair<int,int>> tdc_pairs = {
        {0,1},{2,3},{4,5},{6,7},{8,9},
        {10,11},{12,13},{14,15},{16,17}
    };

    // helper: least-squares tangency solver using normal equations (3x3)
    auto fit_tangent_line = [&](const array<double,6> &xs,
                                const array<double,6> &ys,
                                const array<double,6> &rs,
                                double &out_a, double &out_b, double &out_c) -> bool
    {
        double ATA[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        double ATB[3] = {0,0,0};
        for (int i=0;i<6;++i) {
            double xi = xs[i], yi = ys[i], ri = rs[i];
            ATA[0][0] += xi*xi; ATA[0][1] += xi*yi; ATA[0][2] += xi;
            ATA[1][0] += yi*xi; ATA[1][1] += yi*yi; ATA[1][2] += yi;
            ATA[2][0] += xi;     ATA[2][1] += yi;     ATA[2][2] += 1.0;
            ATB[0] += xi * ri;
            ATB[1] += yi * ri;
            ATB[2] += 1.0 * ri;
        }
        double det = ATA[0][0]*(ATA[1][1]*ATA[2][2]-ATA[1][2]*ATA[2][1])
                   - ATA[0][1]*(ATA[1][0]*ATA[2][2]-ATA[1][2]*ATA[2][0])
                   + ATA[0][2]*(ATA[1][0]*ATA[2][1]-ATA[1][1]*ATA[2][0]);
        if (fabs(det) < 1e-12) return false;
        double M[3][4];
        for (int i=0;i<3;++i) {
            for (int j=0;j<3;++j) M[i][j] = ATA[i][j];
            M[i][3] = ATB[i];
        }
        for (int i=0;i<3;++i) {
            int piv = i;
            for (int r=i+1;r<3;++r) if (fabs(M[r][i]) > fabs(M[piv][i])) piv = r;
            if (fabs(M[piv][i]) < 1e-15) return false;
            if (piv!=i) for (int c=i;c<4;++c) swap(M[i][c], M[piv][c]);
            double div = M[i][i];
            for (int c=i;c<4;++c) M[i][c] /= div;
            for (int r=0;r<3;++r) if (r!=i) {
                double fac = M[r][i];
                for (int c=i;c<4;++c) M[r][c] -= fac * M[i][c];
            }
        }
        double sol[3] = { M[0][3], M[1][3], M[2][3] };
        double a = sol[0], b = sol[1], c = sol[2];
        double norm = sqrt(a*a + b*b);
        if (norm == 0.0) return false;
        out_a = a / norm;
        out_b = b / norm;
        out_c = c / norm;
        return true;
    };

    auto distance_point_line = [&](double a,double b,double c,double x0,double y0)->double{
        return fabs(a*x0 + b*y0 + c);
    };

    // Build index of triggerledge range
    int trigger_min = INT_MAX, trigger_max = INT_MIN;
    for (auto &h: all_hits) {
        trigger_min = min(trigger_min, h.triggerledge);
        trigger_max = max(trigger_max, h.triggerledge);
    }
    if (trigger_min==INT_MAX) {
        cerr << "No hits found\n";
        return 1;
    }
    cout << "Triggerledge range: " << trigger_min << " .. " << trigger_max << "\n";

    vector<SavedHit> out_rows;
    int track_id = 0;

    // windows
    vector<int> windows;
    for (int w=trigger_min; w<=trigger_max; w += WINDOW_SIZE) windows.push_back(w);
    int win_i = 0;
    for (int w0 : windows) {
        ++win_i;
        int w1 = w0 + WINDOW_SIZE;
        cout << "\nProcessing window " << win_i << "/" << windows.size() << ": " << w0 << " - " << w1 << "\n";
        vector<Hit> window_hits;
        window_hits.reserve(4096);
        for (auto &h: all_hits) if (h.triggerledge >= w0 && h.triggerledge < w1) window_hits.push_back(h);
        if (window_hits.empty()) { cout << "  no hits\n"; continue; }

        unordered_map<int, unordered_map<int, vector<const Hit*>>> map_hits;
        map_hits.reserve(32);
        for (auto &h: window_hits) {
            map_hits[h.TDCID][h.CHNLID].push_back(&h);
        }

        // GLOBAL best per top-layer hit across both iterations for this window
        struct BestFit {
            array<const Hit*,6> tube_ptrs;
            array<double,6> xs;
            array<double,6> ys;
            array<double,6> residuals;
            double a,b,c;
            double chi2ndf;
        };
        unordered_map<string, BestFit> global_best_top_map;
        global_best_top_map.reserve(2048);

        // Two iterations
        for (int iteration=1; iteration<=2; ++iteration) {
            array<int,3> layer_offsets;
            array<double,6> signs;
            if (iteration==1) {
                layer_offsets = {0,8,16};
                signs = {+1.0, -1.0, +1.0, +1.0, -1.0, +1.0};
            } else {
                layer_offsets = {0,7,16};
                signs = {-1.0, +1.0, -1.0, -1.0, +1.0, -1.0};
            }

            for (auto &pair_tdcs : tdc_pairs) {
                int t0 = pair_tdcs.first;
                int t1 = pair_tdcs.second;
                if (map_hits.find(t0)==map_hits.end() || map_hits.find(t1)==map_hits.end()) continue;
                auto &chmapA = map_hits[t0];
                auto &chmapB = map_hits[t1];

                for (int base=0; base<8; ++base) {
                    int chA_bot = base + layer_offsets[0];
                    int chA_med = base + layer_offsets[1];
                    int chA_top = base + layer_offsets[2];
                    int chB_bot = chA_bot;
                    int chB_med = chA_med;
                    int chB_top = chA_top;

                    if (chmapA.find(chA_bot)==chmapA.end() || chmapA.find(chA_med)==chmapA.end() || chmapA.find(chA_top)==chmapA.end()) continue;
                    if (chmapB.find(chB_bot)==chmapB.end() || chmapB.find(chB_med)==chmapB.end() || chmapB.find(chB_top)==chmapB.end()) continue;

                    auto &arrA_bot = chmapA[chA_bot];
                    auto &arrA_med = chmapA[chA_med];
                    auto &arrA_top = chmapA[chA_top];
                    auto &arrB_bot = chmapB[chB_bot];
                    auto &arrB_med = chmapB[chB_med];
                    auto &arrB_top = chmapB[chB_top];

                    auto &geoA = geometry_tdcs[t0];
                    auto &geoB = geometry_tdcs[t1];

                    // nested loops (cartesian product)
                    for (const Hit* hA_top : arrA_top) {
                        for (const Hit* hB_top : arrB_top) {
                            for (const Hit* hA_bot : arrA_bot) {
                                for (const Hit* hA_med : arrA_med) {
                                    for (const Hit* hB_bot : arrB_bot) {
                                        for (const Hit* hB_med : arrB_med) {
                                            array<const Hit*,6> tube_ptrs = {hA_bot, hA_med, hA_top, hB_bot, hB_med, hB_top};
                                            array<double,6> xs, ys, rs;
                                            for (int i=0;i<6;++i) {
                                                const Hit* ph = tube_ptrs[i];
                                                int tdc = ph->TDCID;
                                                int ch = ph->CHNLID;
                                                Point p = (tdc==t0) ? geoA[ch] : geoB[ch];
                                                xs[i] = p.first;
                                                ys[i] = p.second;
                                                rs[i] = ph->drift_radius * signs[i];
                                            }
                                            double a,b,c;
                                            bool ok = fit_tangent_line(xs, ys, rs, a,b,c);
                                            if (!ok) continue;
                                            array<double,6> residuals;
                                            double chi2 = 0.0;
                                            for (int i=0;i<6;++i) {
                                                double d = distance_point_line(a,b,c,xs[i],ys[i]);
                                                double res = d - fabs(rs[i]);
                                                residuals[i] = res;
                                                chi2 += res*res;
                                            }
                                            double ndf = 6 - 3;
                                            double chi2ndf = chi2 / ndf;
                                            if (chi2ndf > CHI2NDF_CUT) continue;

                                            // canonical top key: ONLY identify by the top even hit properties (tdc,ch,eventid,triggerledge)
                                            string key = to_string(hA_top->TDCID) + "_" + to_string(hA_top->CHNLID) + "_" +
                                                         to_string(hA_top->eventid) + "_" + to_string(hA_top->triggerledge);

                                            auto it = global_best_top_map.find(key);
                                            if (it == global_best_top_map.end() || chi2ndf < it->second.chi2ndf) {
                                                BestFit bf;
                                                bf.tube_ptrs = tube_ptrs;
                                                bf.xs = xs; bf.ys = ys; bf.residuals = residuals;
                                                bf.a = a; bf.b = b; bf.c = c; bf.chi2ndf = chi2ndf;
                                                global_best_top_map[key] = std::move(bf);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } // end product
                } // end base
            } // end tdc_pairs
        } // end two iterations

        // Save global bests for this window (one entry per top_key)
        int saved = 0;
        for (auto &kv : global_best_top_map) {
            auto &bf = kv.second;
            double tavg = 0.0;
            for (int i=0;i<6;++i) tavg += bf.tube_ptrs[i]->triggerledge;
            tavg /= 6.0;
            for (int i=0;i<6;++i) {
                const Hit* ph = bf.tube_ptrs[i];
                SavedHit sh;
                sh.track_id = track_id;
                sh.TDCID = ph->TDCID;
                sh.CHNLID = ph->CHNLID;
                sh.eventid = ph->eventid;
		sh.drift_time = ph->drift_time;
		sh.corr_time = ph->corr_time;
		sh.adc_time = ph->adc_time;
                sh.triggerledge = ph->triggerledge;
                sh.Dt = ph->triggerledge - tavg;
                sh.x = bf.xs[i];
                sh.y = bf.ys[i];
                sh.drift_radius = ph->drift_radius;
                sh.residual = bf.residuals[i];
                sh.a = bf.a; sh.b = bf.b; sh.c = bf.c;
                sh.chi2ndf = bf.chi2ndf;
                out_rows.push_back(sh);
            }
            cout << "    Saved BEST track " << track_id << " χ2/ndf=" << bf.chi2ndf << "\n";
            ++track_id;
            ++saved;
        }
        cout << "Window saved " << saved << " best tracks\n";
    } // end windows
    
    // ------------------------------------------------------------
// FINAL GLOBAL DEDUPLICATION: ensure only one best track per hA_top
// ------------------------------------------------------------
{
    cout << "\nPerforming global final deduplication across all windows...\n";

    struct TrackBundle {
        vector<SavedHit> hits;   // the 6 hits
        double chi2ndf;
    };

    unordered_map<int, TrackBundle> tracks_by_id;
    tracks_by_id.reserve(out_rows.size()/6);

    // Group hits by track_id
    for (auto &h : out_rows) {
        tracks_by_id[h.track_id].hits.push_back(h);
    }

    // Determine top key per track and find global best
    unordered_map<string, TrackBundle> best_global_top;

    for (auto &kv : tracks_by_id) {
        int tid = kv.first;
        auto &bundle = kv.second;

        if (bundle.hits.size() != 6)
            continue; // should never happen but safety

        // Find the TOP hit (largest y)
        const SavedHit* top_hit = &bundle.hits[0];
        for (auto &h : bundle.hits) {
            if (h.y > top_hit->y) top_hit = &h;
        }

        // Build canonical top key
        string key = to_string(top_hit->TDCID) + "_" +
                     to_string(top_hit->CHNLID) + "_" +
                     to_string(top_hit->eventid) + "_" +
                     to_string(top_hit->triggerledge);

        // Find chi2ndf
        double chi2ndf = bundle.hits[0].chi2ndf;

        // Replace only if better
        auto it = best_global_top.find(key);
        if (it == best_global_top.end() || chi2ndf < it->second.chi2ndf) {
            TrackBundle newb;
            newb.hits = bundle.hits;
            newb.chi2ndf = chi2ndf;
            best_global_top[key] = std::move(newb);
        }
    }

    // Rebuild out_rows
    vector<SavedHit> filtered;
    filtered.reserve(best_global_top.size() * 6);
    for (auto &kv : best_global_top) {
        for (auto &h : kv.second.hits) {
            filtered.push_back(h);
        }
    }

    cout << "Before global dedup: " << out_rows.size()/6 << " tracks\n";
    cout << "After global dedup:  " << filtered.size()/6 << " tracks\n";

    out_rows.swap(filtered);
}

    // write CSV
    cout << "\nSaving output CSV: " << OUTPUT_CSV << "\n";
    ofstream fout(OUTPUT_CSV);
    if (!fout.is_open()) { cerr << "Cannot open output file\n"; return 1; }
    fout << "track_id,TDCID,CHNLID,eventid,drift_time,corr_time,adc_time,triggerledge,Dt,x,y,drift_radius,residual,a,b,c,chi2ndf\n";
    for (auto &r : out_rows) {
        fout << r.track_id << "," << r.TDCID << "," << r.CHNLID << "," << r.eventid << "," << r.drift_time << "," << r.corr_time << "," << r.adc_time << "," << r.triggerledge << ",";
        fout << std::fixed << setprecision(6) << r.Dt << ",";
        fout << std::fixed << setprecision(6) << r.x << "," << r.y << ",";
        fout << std::fixed << setprecision(6) << r.drift_radius << "," << r.residual << ",";
        fout << r.a << "," << r.b << "," << r.c << "," << r.chi2ndf << "\n";
    }
    fout.close();
    cout << "Done. Wrote " << out_rows.size() << " rows (" << (out_rows.size()/6) << " tracks)\n";
    return 0;
}

