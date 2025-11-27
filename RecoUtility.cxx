#include "MuonReco/RecoUtility.h"
#include <iostream>

namespace MuonReco {

  RecoUtility::RecoUtility() {
    MIN_HITS_NUMBER     = 6;
    MAX_HITS_NUMBER     = 12;
    MAX_TIME_DIFFERENCE = 8.0*25.0; // 200 ns is 8 coarse time steps

    MIN_CLUSTER_SIZE    = 2;
    MAX_CLUSTER_SIZE    = 216;

    MIN_CLUSTERS_PER_ML = 1;
    MAX_CLUSTERS_PER_ML = 1;

  }

  RecoUtility::RecoUtility(ParameterSet ps) : RecoUtility() {
    Configure(ps);
  }

  void RecoUtility::Configure(ParameterSet ps) {
    CHECK_TRIGGERS      = ps.getBool  ("CHECK_TRIGGERS",      1,        0);
    IS_PHASE2_DATA      = ps.getBool  ("IS_PHASE2_DATA",      0,        0);
    WIDTHSEL            = ps.getInt   ("WIDTHSEL",            2,        0);
    IS_RELATIVE_DATA    = ps.getBool  ("IS_RELATIVE_DATA",    0,        0);  
    ADC_NOISE_CUT       = ps.getDouble("ADC_NOISE_CUT",       0,        0);
    

    MIN_HITS_NUMBER     = ps.getInt   ("MIN_HITS_NUMBER",     6,        0);
    MAX_HITS_NUMBER     = ps.getInt   ("MAX_HITS_NUMBER",     12,       0);
    MAX_TIME_DIFFERENCE = ps.getDouble("MAX_TIME_DIFFERENCE", 8.0*25.0, 0);

    MIN_CLUSTER_SIZE    = ps.getInt   ("MIN_CLUSTER_SIZE",    2,        0);
    MAX_CLUSTER_SIZE    = ps.getInt   ("MAX_CLUSTER_SIZE",    216,      0);

    MIN_CLUSTERS_PER_ML = ps.getInt   ("MIN_CLUSTERS_PER_ML", 1,        0);
    MAX_CLUSTERS_PER_ML = ps.getInt   ("MAX_CLUSTERS_PER_ML", 1,        0);
    
    TRIGGER_OFFSET      = ps.getDouble("TRIGGER_OFFSET",      0,        0);

    SIG_VOLTAGE_INVERT  = ps.getInt   ("SIG_VOLTAGE_INVERT",  0,        0);
    TRG_VOLTAGE_INVERT  = ps.getInt   ("TRG_VOLTAGE_INVERT",  0,        0);
    std::cout<<"Configure IS_PHASE2_DATA="<<IS_PHASE2_DATA<<std::endl;
    std::cout<<"Configure ADC_NOISE_CUT="<<ADC_NOISE_CUT<<std::endl;
  }

  bool RecoUtility::IsPhase2Data(){return IS_PHASE2_DATA;}

  bool RecoUtility::CheckEvent(Event e, int* status) {
    
    // need precisely one trigger for data
    if (CHECK_TRIGGERS && (IS_PHASE2_DATA==0) && e.TriggerHits().size() != 1) {
      //std::cout << "Rejecting event for not having precisely one trigger (has " << e.TriggerHits().size() << ")" << std::endl;
      *status = 1;
      return NOTPASTEVENTCHECK;
    }
    // need at least 4 wire hits
    if (e.WireHits().size() < MIN_HITS_NUMBER) {
      //std::cout << "Rejecting event for having too few hits (has " << e.WireHits().size() << ")" << std::endl;
      //std::cout << "Event has " << e.WireSignals().size() << " wire signals" << std::endl;
      *status = 2;
      return NOTPASTEVENTCHECK;
    }
    
    if (e.WireHits().size() > MAX_HITS_NUMBER) {
      //std::cout << "Rejecting event for having too many hits (has " << e.WireHits().size() << ")" << std::endl;
      *status = 3;
      return NOTPASTEVENTCHECK;
    }

    int nML0 = 0;
    int nML1 = 0;
    for (Cluster c : e.Clusters()) {
      if (Geometry::MultiLayer(c.Hits().at(0).Layer()) == 0)
        nML0++;
      if (Geometry::MultiLayer(c.Hits().at(0).Layer()) == 1) 
        nML1++;

      // check that no cluster is too large or small
      if (c.Size() < MIN_CLUSTER_SIZE || c.Size() > MAX_CLUSTER_SIZE) {
        //std::cout << "Rejecting event for clusters of wrong size present" << std::endl;
        *status = 4;
        return NOTPASTEVENTCHECK;      
      }
    }

    // check that there is the right number of clusters in each multilayer
    if (nML0 < MIN_CLUSTERS_PER_ML || nML0 > MAX_CLUSTERS_PER_ML ||
      nML1 < MIN_CLUSTERS_PER_ML || nML1 > MAX_CLUSTERS_PER_ML) {
      //std::cout << "Rejecting event for wrong number of clusers in each multilayer" << std::endl;
      //std::cout << "nML0 = " << nML0 << " nML1 = " << nML1 << std::endl;
      //std::cout << "N hits = " << e.WireHits().size() << std::endl;
      *status = 5;
      return NOTPASTEVENTCHECK;
    }

    // need the maximum time difference to be less than max time difference
    double max_time = -1e100;
    double min_time = 1e100;
    double time_diff = 0;

    if (e.WireHits().size()>2){
      for (int i = 0; i < e.WireHits().size(); i++) {
        time_diff = e.WireHits().at(i).DriftTime();
        max_time = max_time < time_diff ? time_diff : max_time;
        min_time = min_time > time_diff ? time_diff : min_time;       
      }
    }

    if ((max_time - min_time <= MAX_TIME_DIFFERENCE) && (max_time - min_time >= 0)) {      
      *status = 0;
      return PASTEVENTCHECK;
    }
    else {
      //std::cout << "Rejecting event for hit timing" << std::endl;
      *status = 6;
      return NOTPASTEVENTCHECK;
    }
  } // end check event


  void RecoUtility::DoHitClustering(Event *e) {
    Bool_t merge = kTRUE;
    std::vector<Cluster> clusters;

    for (auto hit : e->WireHits()) {
      clusters.push_back(Cluster(hit));
    }

    while (merge && clusters.size() >=2) {
      merge = kFALSE;
      for (int i = 0; i < clusters.size(); i++) {
        for (int j = i+1; j < clusters.size(); j++) {
          if (Geometry::AreAdjacent(clusters.at(i),clusters.at(j))) {
            clusters.at(i).Merge(clusters.at(j));
            clusters.erase(clusters.begin()+j);
            merge = kTRUE;
            break;
          }
        } // end for: c2
        if (merge) break;
      } // end for: c1
    } // end while: merge

    // modify event
    for (auto c : clusters) {
      if (c.Size() >= MIN_CLUSTER_SIZE) {
        e->AddCluster(c);
      }
    }
  } //DoHitClustering

  int RecoUtility::DoHitFinding(Event *e, TimeCorrection* tc, Geometry& geo) {
    // std::cout<<"Hit Finding ADC_NOISE_CUT="<<ADC_NOISE_CUT<<std::endl;

    if (IS_PHASE2_DATA){
      return DoHitFindingPhase2(e,   tc, geo, ADC_NOISE_CUT, IS_RELATIVE_DATA);
    }
    else{
      return DoHitFindingPhase1(e,   tc, geo);
    }
  }
    

  int RecoUtility::DoHitFindingPhase1(Event *e, TimeCorrection* tc, Geometry& geo) {
    // Hit h;
    // double adc_time = 0.0;
    // int layer, column;
    // double hx, hy;
    // int nhits = 0;
    // // do trigger hit finding
    // for (auto trig : e->TrigSignals()) {
    //   if (trig.Type() == Signal::RISING) {
    
    // for (auto trig2 : e->TrigSignals()) {
    //   if (trig2.Type() == Signal::FALLING && trig.SameTDCChan(trig2)) {
    //     adc_time = trig2.Time() - trig.Time();
    //     break;
    //   }
    // }
    // if (TRG_VOLTAGE_INVERT) adc_time *= -1.0;
    // // have found adc time for this trigger
    // // construct a hit and push back onto some hit vector
    // if (trig.IsFirstSignal() && adc_time != 0.0 && trig.Time() > (290+TRIGGER_OFFSET) && trig.Time() < (350+TRIGGER_OFFSET)) {
    //   h = Hit(trig.Time(), adc_time, trig.Time(), trig.Time(), trig.TDC(), trig.Channel(), -1, -1, -1, -1);
    //   e->AddTriggerHit(h);
    // }
    //   }
    // } // end for: t_iter
    
    // // do signal hit finding, can only do if triggers exist
    // if (e->TrigSignals().size() != 0) {
    //   Signal selectTrigger = e->TrigSignals().at(0);
    //   double drift_time, corr_time;
    //   Bool_t isFirst = 0;
    //   for (auto sig : e->WireSignals()) {
    //     if (sig.Type() == Signal::RISING) {
    //       drift_time = sig.Time() - selectTrigger.Time();
          
    //       adc_time = 0.0;
    //       for (auto sig2 : e->WireSignals()) {
    //         if ((sig2.Type() == Signal::FALLING) && sig.SameTDCChan(sig2)) {
    //           adc_time = sig2.Time() - sig.Time();
    //           isFirst = (SIG_VOLTAGE_INVERT) ? sig2.IsFirstSignal() : sig.IsFirstSignal();
    //           break;
    //         }
    //       } // end for: s_iter2
    //       if (SIG_VOLTAGE_INVERT) {
    //         adc_time *= -1;
    //         drift_time = drift_time - adc_time;
    //       }
    //       if (isFirst && adc_time > 40) {
    //         corr_time = drift_time - tc->SlewCorrection(adc_time);
            
    //         geo.GetHitLayerColumn(sig.TDC(), sig.Channel(), &layer, &column);
    //         geo.GetHitXY(layer, column, &hx, &hy);
    //         h = Hit(sig.Time(), adc_time, drift_time, corr_time, sig.TDC(), sig.Channel(), layer, column, hx, hy);
    //         e->AddSignalHit(h);
    //         nhits++;
    //       }
    //     }  //if(sig.Type() == Signal::RISING)
    //   } // end for: s_iter
    // } // end if: nonzero number of triggers

    // return nhits;
    return 0;
  } // end method: do hit finding


  int RecoUtility::rollover_bindiff_cal(int a, int b, int rollover){
    int bindiff;
    bindiff = a-b;
    if (bindiff > rollover/2)   bindiff -= rollover;
    else if (bindiff < - (rollover/2)) bindiff += rollover;
    return bindiff;
  }


  int RecoUtility::DoHitFindingPhase2(Event *e, TimeCorrection* tc, Geometry& geo, double adc_cut, int relative) {
    Hit h;
    double adc_time   = 0;
    double tdc_time   = 0;
    double drift_time = 0;
    double corr_time  = 0;
    int layer, column;
    double hx, hy;
    int nhits = 0;
    for (auto sig : e->WireSignals()) {
      adc_time = sig.Width()*BINSIZE*WIDTHSEL;
      if(adc_time>=adc_cut){
        // std::cout<<"ADC cut="<<adc_cut<<std::endl;
        tdc_time = sig.LEdge()*BINSIZE;
        if (relative){
          drift_time = sig.LEdge()*BINSIZE;
        }
        else {
          drift_time = rollover_bindiff_cal(sig.LEdge(),e->TrigSignals().at(0).TriggerLEdge(),ROLLOVER)*BINSIZE;
	  //std::cout << "drift_time " << drift_time << std::endl;
	  //std::cout << "ledge" << sig.LEdge() << std::endl;
	  //std::cout << "trigger" << e->TrigSignals().at(0).TriggerLEdge() << std::endl;
          
        }
        corr_time = drift_time - tc->SlewCorrection(adc_time);
        geo.GetHitLayerColumn(sig.TDC(), sig.Channel(), &layer, &column);
        geo.GetHitXY(layer, column, &hx, &hy);
        h = Hit(sig.Time(), adc_time, drift_time, corr_time, sig.TDC(), sig.Channel(), layer, column, hx, hy);
        e->AddSignalHit(h);
        nhits++;

	// === DEBUG: dump event_id and TDC to CSV in chunks of 1m hits ===
        static std::ofstream debugCSV;
        static long long debugHitCounter = 0;
        static int debugFileIndex = 0;

        if (!debugCSV.is_open()) {
          std::string filename = "hits_" + std::to_string(debugFileIndex) + ".csv";
          debugCSV.open(filename);
          debugCSV << "eventid,eventid_t,eventid_ext,CSMID,TDCID,CHNLID,WIDTH,TDC_EVENTID,TDC_BCID,ledge,mode,triggerledge,hdrtrlid,adc_time, drift_time, corr_time,layer, column, hx, hy\n";
        }
        
	debugCSV << static_cast<long>(sig.HeaderEID()) << "," << static_cast<long>(sig.TrailerEID()) << "," << static_cast<long>(sig.HeaderEIDext()) <<"," << static_cast<int>(sig.CSMID()) << "," << static_cast<int>(sig.TDC()) << "," << static_cast<int>(sig.Channel()) <<"," << static_cast<int>(sig.Width()) <<"," << static_cast<int>(sig.TDCHeaderEID()) <<"," << static_cast<int>(sig.TDCHeaderBCID()) <<"," << static_cast<long long>(sig.LEdge()) << ","<< static_cast<int>(sig.Mode()) << "," << static_cast<long long>(sig.TriggerLEdge()) <<","<< static_cast<int>(sig.TDCHdrTrlrID()) << "," << adc_time << "," << drift_time <<"," << corr_time <<"," << layer << "," << column << "," << hx <<","<< hy << "\n";
        debugHitCounter++;

        if (debugHitCounter >= 1000000) {
         debugCSV.close();
         debugHitCounter = 0;
         debugFileIndex++;
         std::string filename = "hits_" + std::to_string(debugFileIndex) + ".csv";
         debugCSV.open(filename);
         debugCSV << "eventid,eventid_t,eventid_ext,CSMID,TDCID,CHNLID,WIDTH,TDC_EVENTID,TDC_BCID,ledge,mode,triggerledge,hdrtrlid,adc_time, drift_time, corr_time,layer, column, hx, hy\n";
        }
      } //if(adc_time>=adc_cut)
    } //for (auto sig : e->Signals())
    return nhits;
  }
} 
