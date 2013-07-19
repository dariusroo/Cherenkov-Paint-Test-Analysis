/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Code has been modified for use by SNO+ Collaboration                      *
 * Original Copyright notice below.                                          *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// C Standard Libraries
#include <iostream>
#include <string>
#include <fstream>

// Root Libraries
#include <TFile.h>
#include <TMath.h>
#include <TRandom2.h>
#include <TMathBase.h>
#include <TTree.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TPaveStats.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TLine.h>
#include <TGaxis.h>

// HDF5 Library
#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

    using namespace std;

/* Global Variables */

const int    RANK_OUT = 2;             // Data Rank
const double LOW_LIM = 0.5;            // Upper Limit of Low-Charge Region (in picoCoulomb or "pC")
const int NUM_TEST_TRACES = 100000;    // Number of Traces to Use for Acceptable Window Range Determination

/* Struct: DataCluster
 * * * * * * * * * * *
 * Stores all parameters necessary to extract a trace from an initialized DataSet.
 * Contains data buffer.
 */

struct DataCluster{
  DataSet *dataset;                    // Dataset pointer
  DataSpace dataspace;                 // DataSet's DataSpace
  DataSpace memspace;                  // MemSpace Object for Data Extraction
  hsize_t      offset[RANK_OUT];       // Data Extraction Parameters...
  hsize_t      count[RANK_OUT];
  hsize_t      offset_out[RANK_OUT];
  hsize_t      count_out[RANK_OUT];
  unsigned long trace_length;          // Length of a Scope Trace (each subarray is a trace in sequence mode) ... unsigned long is probably unecessary
  unsigned long n_traces;              // Number of traces in DataSet
  char * data_out;                     // Pointer to Data Buffer
};

typedef struct DataCluster DataCluster;

/* DataCluster Methods */

DataCluster * Init_Data(DataSet *dataset);
int Read_Trace(DataCluster *datacluster, unsigned long trace_index);

/* Pedestal Correction Methods */

//float variancer(float variance, float voltage, unsigned long i);


/* Method: main(int argc, char* argv[])
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * ./ReadData_C <input_file_list.txt> <root_file_output.root> <window_width> <maximum_range>(opt.)
 *
 * <input_file_list.txt> : Text file containing .hdf5 files to be analyzed, spaced by newlines
 * <root_file_output.root> : Histogram output root file (should not exist already)
 * <window_width> : Pedestal correction window in units of time bins
 * <maximum_range> : Maximum pedestal window range to be used in the analysis (bypasses SetRange()) in units of voltage bins
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * Reads in a list of .hdf5 LeCroy WaveRunner 6 Zi Oscilloscope data files and outputs histograms of various quantities. 
 * Implements an event-specific pedestal correction.
 */

int main (int argc, char* argv[])
{
    int window_width = atoi(argv[3]); // Width of pedestal correction region
 
   try
     {

     // Attribute Variables
     Attribute horiz_interval;
     Attribute vertical_gain;
     Attribute horiz_offset;
     Attribute vertical_offset;
     Attribute max_value;
     Attribute min_value;
     
     // Attribute Value Variables
     double dx,dy,xoffset,yoffset;
     double Vmin, Vmax;

     // Analysis Variables 
     TH1F *charges = new TH1F("LED Charge","",5000,-150,200);    // Charge Spectrum Post-Trigger 
     TH1F *ncharges = new TH1F("Noise Charge","",5000,-150,200); // Charge Spectrum Pre-Trigger
     TH1F *pedestals = new TH1F("LED Pedestal Correction","",1000,-0.01,0.01);  // Pedestal Correction for LED
     TH1F *npedestals = new TH1F("Noise Pedestal Correction","",1000,-0.01,0.01);  // Pedestal Correction for noise
     TH1F *cpedestals = new TH1F("Corrected LED Pedestal Correction","",1000,-0.01,0.01); // Corrected Pedestal Correction for LED
     TH1F *cnpedestals = new TH1F("Corrected Noise Pedestal Correction","",1000,-0.01,0.01); // Corrected Pedestal Correction for LED
     TH1F *vtrace = new TH1F("Lower Variance Traces","",10002,0,1);            // Traces with Variance between 0.0003 and 0.0004 
     TH1F *htrace = new TH1F("Higher Variance Traces","",10002,0,1);            // Traces with Variance between above 0.00014 
     TH1F *gped =  new TH1F("LED Pedestal Graph","",1000000,0,1);    // Pedestal Plot For Time-Evolution
     TH1F *gcharge =  new TH1F("Charge Graph","",1000000,0,1);   // Charge Plot For Time-Evolution
     TH1F *voltages = new TH1F("Corrected Voltages","",30000,-0.5,0.1);  // Pedestal-Corrected Voltages (mV)
     TH1F *voltages1 = new TH1F("Corrected Voltages1","",30000,-0.5,0.1);  // Pedestal-Corrected Voltages (mV)
     TH1F *variances = new TH1F("Variance Charge","",2000000,-0.5,0.5); // Variance 
     TH1F *nvariances = new TH1F("Variance NCharge","",2000000,-0.5,0.5); // Variance 

     unsigned long i;            // One size fits all
     unsigned long p;
     unsigned long k;
     unsigned long b; 
     int f = 0.0;
     int q = 0.0; 
     int g = 0.0; 

     int tcount = 0;             // Trace count

     int lower_count=0.0; 
     int higher_count=0.0; 
 

      // Overriding SetRange()
     // if (argc == 5){ 
     //	maximum_range = atoi(argv[4]);
     //	cout << endl << "Acceptable Pedestal Window Range Set @ "<< maximum_range * dy  <<" mV for this analysis."<< endl; 
     //  }
      
      for(p=0; p < 1; p++){

	// Data Reading Variables
	string filename;
	H5File file;
	DataSet dataset;

	ifstream ifs ( argv[1] , ifstream::in ); // Open File List Stream
	ifs >> filename;

	while (ifs.good()){  // While Stream Is Open, Analyze New Files.
      
	  cout<<filename<<endl;      
	  file.openFile(filename, H5F_ACC_RDONLY ); // Open HDF5 File
	  ifs >> filename;

	  dataset = file.openDataSet("channel1" );  // Open HDF5 DataSet
	
	  // Reading in Attributes

	  horiz_interval  = dataset.openAttribute("horiz_interval");
	  vertical_gain  = dataset.openAttribute("vertical_gain");
	  horiz_offset  = dataset.openAttribute("horiz_offset");
	  vertical_offset  = dataset.openAttribute("vertical_offset");
	  max_value  = dataset.openAttribute("max_value");
	  min_value  = dataset.openAttribute("min_value");
	
	  horiz_interval.read(PredType::NATIVE_DOUBLE, &dx);
	  vertical_gain.read(PredType::NATIVE_DOUBLE, &dy);
	  horiz_offset.read(PredType::NATIVE_DOUBLE, &xoffset);
	  vertical_offset.read(PredType::NATIVE_DOUBLE, &yoffset);
	  max_value.read(PredType::NATIVE_DOUBLE, &Vmax);
	  min_value.read(PredType::NATIVE_DOUBLE, &Vmin);
	
	  // Processing Attributes

	  int window_offset = -(int)(xoffset/dx); 
	  cout << "Trigger Offset: " << window_offset*dx << " seconds (" << window_offset << " points)." << endl;
	 

	 

	  //cout << "Minimum voltage: "<< ((float)Vmin*dy-yoffset) <<" mV, Maximum: "<< ((float)Vmax*dy-yoffset) <<" with offset: "<<yoffset <<  " Vertical gain: " << dy << endl;
	  DataCluster * datacluster = Init_Data(&dataset); // Initialize that DataCluster mess



 
	  // Analysis Variables
      
	  double charge=0.0;
	  float pedestal = 0.0; 
	  float cpedestal =0.0; 
	  float voltage;
	  float voltager; 
	  float variance;  
	  unsigned long integrate_length = TMath::Min((unsigned long)(window_offset - 2*window_width),datacluster->trace_length - window_offset);
	  cout<<" The window offest is "<<  window_offset << " The window width is "<<window_width<<" The integrate length is "<<integrate_length<<" dy is "<< dy <<" dx is "<< dx <<" x offset is "<<xoffset<< endl; 


	  // For each trace
	  for (unsigned long j = 0; j < datacluster->n_traces; j++)
	    {
	      cout<<"  Analyzing trace "<<j+1<<" out of "
		  <<datacluster->n_traces<<" total traces"<<'\r';
	      cout.flush();

	      Read_Trace(datacluster,j);
	      tcount++;	
      
		// Trace-Specific Pedestal Correction 
		pedestal = TMath::Mean (window_width, &datacluster->data_out[window_offset-window_width])*dy;
		gped->SetBinContent(tcount + 1, pedestal);

		charge=0.0;
		variance=0.0; 
		f = 0.0; 
	

		// Look past the trigger offset


		for(i = window_offset; i < window_offset + window_width ; i++){
		  voltage = ((float)datacluster->data_out[i]*dy-pedestal); //mV
		  voltages->Fill(voltage);
		  variance = variance + voltage*voltage;
		  f++; 
		}

		if (variance < 0.000018){
		  cpedestal =  TMath::Mean (window_width, &datacluster->data_out[window_offset-window_width])*dy;
		}
	      
		for(b = window_offset; b < window_offset + integrate_length ; b++){
		  voltager = ((float)datacluster->data_out[b]*dy-cpedestal); //mV
		  charge=charge+(voltager*((-1000.0*dx*1e9)/75.0));
		}
	    	 
		if (variance > 0.00003 && variance < 0.00004){
		  for(i = 0; i < datacluster->trace_length; i++){                                            
		    vtrace->SetBinContent(i+1,(float)datacluster->data_out[i]*dy-pedestal); 
		  }
		  lower_count++;
		}

		if (variance > 0.00016){
		  for(i = 0; i < datacluster->trace_length; i++){                                            
		    htrace->SetBinContent(i+1,(float)datacluster->data_out[i]*dy-pedestal); 
		  }
		  higher_count++;
		}
	
		// Store distributions
		variances->Fill(variance);
		charges->Fill(charge);
		pedestals->Fill(pedestal);
		cpedestals->Fill(cpedestal); 
		gcharge->SetBinContent(tcount + 1, charge);


	    }
	  file.close();
	} 
      }

      cout << "Number of traces with variance between 0.00003 and 0.00004 (mV^2) is " << lower_count << endl; 
      cout << "Number of traces with variance over 0.00016 (mV^2) is " << higher_count << endl; 
      cout << "Number of Bins for each traces in the post-trigger pedestal window is " << f << endl; 
      //   cout << "Number of Traces in the LED pedestal window that are accepted is " << window_count << endl; 

      for(k=0; k < 1; k++){


	// Data Reading Variables
	string filename;
	H5File file;
	DataSet dataset;

	ifstream ifs ( argv[1] , ifstream::in ); // Open File List Stream
	ifs >> filename;

	while (ifs.good()){  // While Stream Is Open, Analyze New Files.
      
	  cout<<filename<<endl;      
	  file.openFile(filename, H5F_ACC_RDONLY ); // Open HDF5 File
	  ifs >> filename;

	  dataset = file.openDataSet("channel1" );  // Open HDF5 DataSet
	
	  // Reading in Attributes

	  horiz_interval  = dataset.openAttribute("horiz_interval");
	  vertical_gain  = dataset.openAttribute("vertical_gain");
	  horiz_offset  = dataset.openAttribute("horiz_offset");
	  vertical_offset  = dataset.openAttribute("vertical_offset");
	  max_value  = dataset.openAttribute("max_value");
	  min_value  = dataset.openAttribute("min_value");
	
	  horiz_interval.read(PredType::NATIVE_DOUBLE, &dx);
	  vertical_gain.read(PredType::NATIVE_DOUBLE, &dy);
	  horiz_offset.read(PredType::NATIVE_DOUBLE, &xoffset);
	  vertical_offset.read(PredType::NATIVE_DOUBLE, &yoffset);
	  max_value.read(PredType::NATIVE_DOUBLE, &Vmax);
	  min_value.read(PredType::NATIVE_DOUBLE, &Vmin);
	
	  // Processing Attributes

	  int window_offset = -(int)(xoffset/dx); 
	  cout << "Trigger Offset: " << window_offset*dx << " seconds (" << window_offset << " points)." << endl;
	  cout << window_offset << window_width << endl; 



	  //cout << "Minimum voltage: "<< ((float)Vmin*dy-yoffset) <<" mV, Maximum: "<< ((float)Vmax*dy-yoffset) <<" with offset: "<<yoffset <<  " Vertical gain: " << dy << endl;
	  DataCluster * datacluster = Init_Data(&dataset); // Initialize that DataCluster mess



	  // Analysis Variables
      
	  double ncharge = 0.0;
	  float npedestal = 0.0;
	  float cnpedestal = 0.0; 
	  float voltage1;
	  float voltager1; 
	  float nvariance; 
	  unsigned long integrate_length = TMath::Min((unsigned long)(window_offset - 2*window_width),datacluster->trace_length - window_offset);
 

	  // For each trace
	  for (unsigned long j = 0; j < datacluster->n_traces; j++)
	    {
	      cout<<"  Analyzing trace "<<j+1<<" out of "
		  <<datacluster->n_traces<<" total traces"<<'\r';
	      cout.flush();

	      Read_Trace(datacluster,j);
	      tcount++;
 
	   
		// Trace-Specific Pedestal Correction 
		npedestal = TMath::Mean (window_width, datacluster->data_out)*dy;

		nvariance=0.0;
		ncharge=0.0;
		g = 0.0; 
  		// Look before pedestal window in noise region - This and the former should be the same length 
		

		for(i = 0.0; i < window_width; i++){
		  voltage1 = ((float)datacluster->data_out[i]*dy-npedestal); //mV
		  voltages1->Fill(voltage1);
		  nvariance = nvariance + voltage1*voltage1; 
		  g++;
		}
	
		if (nvariance < 0.000018){
		 cnpedestal = TMath::Mean (window_width, datacluster->data_out)*dy; 
		}

		for(b = window_width; b < window_width + integrate_length; b++){
		  voltager1 = ((float)datacluster->data_out[b]*dy-cnpedestal); //mV	  
		  ncharge=ncharge+(voltager1*((-1000.0*dx*1e9)/75.0)); //pC
		}
		// Store distributions
		nvariances->Fill(nvariance);
		npedestals->Fill(npedestal);
		ncharges->Fill(ncharge);
		cnpedestals->Fill(cnpedestal);


		// }
		//	 else { // If window was bad, reject the trace (you could also loop over the noise region to find a better one, but you'd have to think about what that means...)
	      //	reject_count++;
		   //	 }
	}
	  file.close();
	  delete[] datacluster->data_out ;
	  delete[] datacluster ;
      }
      ifs.close();
     }
//  cout << endl;
//   cout << "Rejected " << reject_count << " traces."<< endl;
      cout << "Number of Bins for each traces in the pre-trigger pedestal window is " << g << endl; 
      //  cout << "Number of Traces in the noise pedestal window that are accepeted is " << window_counter << endl; 
      cout << q << endl; 

    // Output Histograms to File
      if (argc > 2){
	TFile f(argv[2],"new");
	charges->Write();
	ncharges->Write();
	pedestals->Write();
	npedestals->Write();
	cnpedestals->Write();
	voltages->Write();
	gped->Write();
	gcharge->Write();
	vtrace->Write();
	htrace->Write();
	voltages1->Write();
	variances->Write();
	nvariances->Write();
      
      }
      //Clean-up time!
       delete charges ;
       delete ncharges ;
       delete pedestals ;
       delete npedestals ;
       delete cnpedestals ; 
       delete voltages ;
       delete gped ;
       delete gcharge ;
       delete vtrace ;
       delete htrace ; 
       delete voltages1 ;
       delete variances ; 
       delete nvariances ; 

  
     }  // end of try block



   // catch failure caused by the H5File operations
   catch( FileIException error )
     {
      error.printError();
      return -1;
     }

   // catch failure caused by the DataSet operations
   catch( DataSetIException error )
   {
      error.printError();
      return -1;
   }

   // catch failure caused by the DataSpace operations
   catch( DataSpaceIException error )
   {
      error.printError();
      return -1;
   }

   // catch failure caused by the DataSpace operations
   catch( DataTypeIException error )
   {
      error.printError();
      return -1;
   }

   return 0;  // successfully terminated
}

const int ABS_RAT = 2; 

/* Method: SetRange(int n_test_traces, DataCluster *datacluster, int max_test_range, int window_width, int window_offset)
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * Determines the maximum pedestal window range by increasing it until 1/ABS_RAT or more of the traces are accepted.
 */

//float variancer(float variance, float voltage, unsigned long i,){
// voltage = ((float)datacluster->data_out[i]*dy-pedestal); 
// variance = variance + voltage*voltage; 
// return (variance,voltage,i); 
//}
/* Method: Init_Data(DataSet *dataset)
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * Returns a DataCluster... 'nuff said, hopefully...
 */

DataCluster * Init_Data(DataSet *dataset){

  DataCluster * datacluster = new DataCluster[1];
  
  /* 	
   * Get dataspace of the dataset.
   */
  datacluster->dataset = dataset;
  datacluster->dataspace = datacluster->dataset->getSpace();
  
  /*
   * Get the dimension size of each dimension in the dataspace and
   * display them.
   */
  
  hsize_t dims_out[2];
  datacluster->dataspace.getSimpleExtentDims( dims_out, NULL);
  datacluster->trace_length = (unsigned long)(dims_out[1]);
  datacluster->n_traces = (unsigned long)(dims_out[0]);
      
  // cout << "Reading " <<  datacluster->n_traces << " traces of length " << datacluster->trace_length << "..." << endl;
  
  // Data Buffer
  datacluster->data_out = new char[datacluster->trace_length]; // Scope data is size char
  for (unsigned long i = 0; i < datacluster->trace_length; i++) datacluster->data_out[i]= 0;

  /*
   * Define hyperslab in the dataset.
   */

  datacluster->offset[0] = 0;
  datacluster->offset[1] = 0;
  datacluster->count[0]  = 1;
  datacluster->count[1]  = datacluster->trace_length;
  datacluster->dataspace.selectHyperslab( H5S_SELECT_SET, datacluster->count, datacluster->offset );
  
  /*
   * Define the memory dataspace.
   */

  hsize_t     dimsm[2];              /* memory space dimensions */
  dimsm[0] = dims_out[0];
  dimsm[1] = dims_out[1];
  datacluster->memspace = DataSpace( RANK_OUT, dimsm );

  /*
   * Define memory hyperslab.
   */

  datacluster->offset_out[0] = 0;
  datacluster->offset_out[1] = 0;
  datacluster->count_out[0]  = 1;
  datacluster->count_out[1]  = datacluster->trace_length;
  datacluster->memspace.selectHyperslab( H5S_SELECT_SET, datacluster->count_out, datacluster->offset_out );
  
  return datacluster;
}

/* Method: Read_Trace(DataCluster *datacluster, unsigned long trace_index)
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * Updates a DataCluster datacluster so that its buffer contains trace number trace_index
 */

int Read_Trace(DataCluster *datacluster, unsigned long trace_index){
  datacluster->offset[0]= (hsize_t)trace_index;
  datacluster->dataspace.selectHyperslab( H5S_SELECT_SET, datacluster->count, datacluster->offset );
  datacluster->memspace.selectHyperslab( H5S_SELECT_SET, datacluster->count_out, datacluster->offset_out );
  datacluster->dataset->read( datacluster->data_out, PredType::NATIVE_CHAR, datacluster->memspace, datacluster->dataspace );
  return 0; // No protection...
}
