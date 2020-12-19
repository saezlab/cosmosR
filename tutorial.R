library(cosmos)

#In order to adapt options to users specification we can load them into a variable 
#that will then be passed to preprocess_COSMOS_signaling_to_metabolism CARNIVAL_options parameter
my_options <- default_CARNIVAL_options()

#Here the user should provide a path to its CPLEX executable (only cplex at the moment, other solvers will be documented soon !)
my_options$solverPath <- "~/Documents/cplex"

#### FORWARD run of COSMOS, to connect signaling to metabolism

#The signaling inputs are the result of footprint based TF and kinase activity estiamtion
#For more info on TF activity estiamtion from transcriptomic data, see:https://github.com/saezlab/transcriptutorial (Especially chapter 4)

test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = toy_sif,
                                                      signaling_data = toy_signaling_input_carnival_vec,
                                                      metabolic_data = toy_metab_input_carnival_vec,
                                                      diff_expression_data = toy_RNA,
                                                      maximum_network_depth = 15,
                                                      remove_unexpressed_nodes = T,
                                                      CARNIVAL_options = my_options
                                                      
)

test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                      CARNIVAL_options = my_options)

metab_to_pubchem_vec <- metab_to_pubchem$name
names(metab_to_pubchem_vec) <- metab_to_pubchem$pubchem

test_result_for <- format_COSMOS_res(test_result_for,
                                     metab_mapping = metab_to_pubchem_vec,
                                     measured_nodes = unique(c(names(toy_metab_input_carnival_vec),names(toy_signaling_input_carnival_vec))),
                                     omnipath_ptm = omnipath_ptm)

View(test_result_for[[1]]) #SIF
View(test_result_for[[2]]) #ATTRIBUTES

#### BACKWARD run of COSMOS, to connect metabolism to signaling

test_back <- preprocess_COSMOS_metabolism_to_signaling(meta_network = toy_sif,
                                                       signaling_data = toy_signaling_input_carnival_vec,
                                                       metabolic_data = toy_metab_input_carnival_vec,
                                                       diff_expression_data = toy_RNA,
                                                       maximum_network_depth = 15,
                                                       remove_unexpressed_nodes = F,
                                                       CARNIVAL_options = my_options
                                                       
)

test_result_back <- run_COSMOS_metabolism_to_signaling(data = test_back,
                                                       CARNIVAL_options = my_options)


test_result_back <- format_COSMOS_res(test_result_back,
                                      metab_mapping = metab_to_pubchem_vec,
                                      measured_nodes = unique(c(names(toy_metab_input_carnival_vec),names(toy_signaling_input_carnival_vec))),
                                      omnipath_ptm = omnipath_ptm)

View(test_result_back[[1]]) #SIF
View(test_result_back[[2]]) #ATTRIBUTES

###Merge forward and backward networks

full_sif <- as.data.frame(rbind(test_result_for[[1]], test_result_back[[1]]))
full_attributes <- as.data.frame(rbind(test_result_for[[2]], test_result_back[[2]]))

full_sif <- unique(full_sif)
full_attributes <- unique(full_attributes)

###Visualise network

#you can give here a signle nodes name or a vector or names. 
#This functions will allow to display the neighboorhood of given nodes

network_plot <- display_node_neighboorhood(central_node = 'PRKACA', sif = full_sif, att = full_attributes, n = 5)

network_plot
