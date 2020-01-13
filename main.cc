#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <fstream>


///////////// PARAMETERS CLASS ///////////////////////////
//

class cParams{

public:
  int grid_size; // Grid size
  int num_tasks; // Number of tasks including task zero
  int comp_cap; // Computational capacity of each org (genome size)

  double mut_rate; // Mutation rate
  int run_time; // Runtime

  double export_eff; //Export efficacy of sharing

  double task_info_a; // Slope of information vs task number linear relationship
  double task_info_b; // Intercept of information vs task number linear relationship

  // Extra parameters for utility functions
  double extr1;
  double extr2;

  // UF ID
  int uf_id;

  //
  int print_every;

  cParams(){ // Constructor initializes all values using user input

		std::ifstream inFile;
		inFile.open("params.cfg");
		if (!inFile) {
			std::cerr << "Unable to open parameters file params.cfg\n";
			exit(1);   // call system to stop
		}
    std::string param_name;
    double i;

		while(inFile >> param_name >> i){
    	if(param_name=="GRID_SIZE"){
        std::cout<<"Grid size found!" + std::to_string(i) +"\n";
        grid_size=(int) i;
      }
    	if(param_name=="NUM_TASKS"){
        std::cout<<"Number of tasks found!"+ std::to_string(i) +"\n";
        num_tasks = (int) i;
      }
    	if(param_name=="COMP_CAP"){
        std::cout<<"Computational capacity found!"+ std::to_string(i) +"\n";
        comp_cap = (int)i;
      }
    	if(param_name=="EXP_EFF"){
        std::cout<<"Export efficiency found!"+ std::to_string(i) +"\n";
        export_eff = i;
      }
    	if(param_name=="MUT_RATE"){
        std::cout<<"Mutation Rate found!"+ std::to_string(i) +"\n";
        mut_rate = i;
      }
    	if(param_name=="RUN_TIME"){
        std::cout<<"Run Time found!"+ std::to_string(i) +"\n";
        run_time = (int) i;
      }
    	if(param_name=="TASK_INFO_A"){
        std::cout<<"Task to Info A found!"+ std::to_string(i) +"\n";
        task_info_a = i;
      }
    	if(param_name=="TASK_INFO_B"){
        std::cout<<"Task to Info B found!"+ std::to_string(i) +"\n";
        task_info_b = i;
      }
    	if(param_name=="UF_ID"){
        std::cout<<"UF ID found!"+ std::to_string(i) +"\n";
        uf_id = (int) i;
      }
    	if(param_name=="EXTR1"){
        std::cout<<"Extra Parameter 1 found!"+ std::to_string(i) +"\n";
        extr1 = i;
      }
    	if(param_name=="EXTR2"){
        std::cout<<"Extra Parameter 2 found!"+ std::to_string(i) +"\n";
        extr2 = i;
      }
			if(param_name=="PRINT_EVERY"){
				std::cout<<"Print every update parameter found!"+ std::to_string(i) +"\n";
				print_every = i;
			}
    };
  }
  double TaskNumToInfo(int task_no){
  	if (task_no != 0){
      return task_info_a * (double)task_no + task_info_b;
  	}
  	else {
  		return 0.0;
  	}
  }
};


/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// UTILITY FUNCTIONS //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

double UF(int id,std::vector<double> info,double extr1, double extr2){
  if(id==0){
    double fn0 = 0.0;
    for ( int i = 0; i < info.size(); i = i + 1){
      fn0 = fn0 + info[i];
    }
    return fn0;
  }

  if(id==1){
    double fn1 = 0.0;
    double tobeadded = 0.0;
    for ( int i = 0; i < info.size(); i = i + 1){
      tobeadded = info[i];
      if (info[i] > 0) {
        for ( int j = 0; j < i; j = j + 1){
          if (info[j] > 0) {
            tobeadded = tobeadded + extr1;
          }
        }
      }
       fn1 = fn1 + tobeadded;
    }
    return fn1;
  }

  if(id==2){
    double fn2 = 0.0;
    double tobeadded = 0.0;
    for ( int i = 0; i < info.size(); i++){
      tobeadded = info[i];
       if (info[i] > 0) {
        for ( int j = 0; j < i; j = j + 1){
          tobeadded = tobeadded + extr1*info[j];
        }
       }
      fn2 = fn2 + tobeadded;
    }
    return fn2;
  }

  if(id==3){
    double fn3 = 0.0;
    double tobeadded = 0.0;
    double net_insight;
    for ( int i = 0; i < info.size(); i = i + 1){
      tobeadded = info[i];
      if (info[i] > 0) {
        for ( int j = 0; j < i; j = j + 1){
          if (info[j] > 0) {
            for ( int k = 0; k <= j; k = k + 1){
              net_insight = pow(extr1, i -k);
              tobeadded = tobeadded + net_insight*info[k];
            }
          }
        }
       }
      fn3 = fn3 + tobeadded;
    }
    return fn3;
  }
}
};


///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// RUNTIME ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
int main ()
{
  srand(static_cast<unsigned int>(time(nullptr))); // SEED as SYSTEM CLOCK
  // Get all required parameters from user

	std::ofstream outputFile;
	outputFile.open("output.dat");

  std::cout << "MAIN WORLD PARAMETERS\n";

  cParams params; // All parameter values store in run_params now

  // Initialize grid and grid cells.
  cGrid grid(params.grid_size,params.comp_cap,params.num_tasks);


  ////// RUN!
  for(int cur_upd = 1; cur_upd<params.run_time ; cur_upd++){
		std::string pix = std::to_string(cur_upd);
    std::cout << "The current run is " + pix + "\n";


    /////Iterate through all cells and calculate info and rec std::vectors
    for(int cell_idx=0; cell_idx < params.grid_size*params.grid_size; cell_idx++){

      cGridCell* cur_cell = grid.GetCell(cell_idx);

      // For this cell, execute positive genomic instructions to form initial info std::vector and get the information it has to send
      std::vector<double> to_dist = cur_cell->ExecuteGenomeTasks(params);

      // distribute to_dist to eight neighbors

      grid.SendToNeighbors(cell_idx,to_dist,params);
		}


	  /////// Add Received to info and calculate fitnesses
	  for(int cell_idx=0; cell_idx < params.grid_size*params.grid_size; cell_idx++){
	    cGridCell* cur_cell = grid.GetCell(cell_idx);
	    cur_cell->AddReceivedToInfo();
	    cur_cell->CalculateFitness(params); /// CHANGE 0 to UF ID
	  }

	  ///////// PRINTING LOOP
	if(cur_upd%params.print_every==0){
		outputFile << "\nCURRENT_TIME "+std::to_string(cur_upd)+"\n";
		outputFile << "The entropy is "+std::to_string(Entropy(grid, params))+"\n";
		outputFile << "cell_id     genome    fitness    EI-index\n";
	  for(int cell_idx=0; cell_idx < params.grid_size*params.grid_size; cell_idx++){
			std::string to_write = ""+std::to_string(cell_idx)+" ";
		cGridCell* cur_cell = grid.GetCell(cell_idx);
		std::vector<int> cur_genome = cur_cell->GetGenome();
		std::vector<int>::iterator it;
		for(it = cur_genome.begin(); it != cur_genome.end(); it++) {
		  to_write += "("+std::to_string(*it)+")";
		}
		double fitness_x = cur_cell->GetFitness();
		to_write += " "+std::to_string(fitness_x)+" ";
			to_write += std::to_string(CalcEIRatio(cur_genome))+"\n";
			outputFile << to_write;
	}
	}

    ////// Reproduce!
    cGrid new_grid(params.grid_size,params.comp_cap,params.num_tasks);

    for(int cell_idx=0; cell_idx < params.grid_size*params.grid_size; cell_idx++){

      std::vector<int> genome_at_idx = grid.ReproGenome(cell_idx,params);
      new_grid.GetCell(cell_idx)->SetGenome(genome_at_idx);
    }

    for(int cell_idx=0; cell_idx < params.grid_size*params.grid_size; cell_idx++){
      grid.GetCell(cell_idx)->SetGenome(new_grid.GetCell(cell_idx)->GetGenome());
    }

    ////// Clear info std::vectors, rec std::vectors and fitnesses.
    for(int cell_idx=0; cell_idx < params.grid_size*params.grid_size; cell_idx++){
      cGridCell* cur_cell = grid.GetCell(cell_idx);
      cur_cell->SetFitness(0);
      cur_cell->ClearReceived();
      cur_cell->ClearInfo();
    }
	}
	outputFile.close();
  return 0;
}
