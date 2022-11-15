#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

int num_of_CT1_Ligands;
int num_of_CT1_Receptors;
int num_of_CT2_Ligands;
int num_of_CT2_Receptors;
int num_of_CT1_CMs;
int num_of_CT2_CMs;
int num_of_LR_pairs;

int **interaction_array;

typedef struct Expression{
    char original_name[20];
    char gene_ID[20];
    int *CM_exp;
}Expression;

// ligand-recrptor pair from database
typedef struct Pair{
    char pair_ori_name[40];
    char pair_name[40];
    char cluster_ligand[20];
    char cluster_receptor[20];
    int cluster_ID;
    int hit;
}Pair;

// cluster name or ID
typedef struct cluster_ID{
    char Col_Name[20];
}clster_ID;

Expression *CT1_Cond1_Ligand_exp_table;
Expression *CT1_Cond1_Receptor_exp_table;
Expression *CT1_Cond2_Ligand_exp_table;
Expression *CT1_Cond2_Receptor_exp_table;
Expression *CT2_Cond1_Ligand_exp_table;
Expression *CT2_Cond1_Receptor_exp_table;
Expression *CT2_Cond2_Ligand_exp_table;
Expression *CT2_Cond2_Receptor_exp_table;

Pair *LR_table;
Pair *Cond1_CT2_ligand_to_CT1_receptor_pair_list;
Pair *Cond2_CT2_ligand_to_CT1_receptor_pair_list;
Pair *Cond1_CT1_ligand_to_CT2_receptor_pair_list;
Pair *Cond2_CT1_ligand_to_CT2_receptor_pair_list;

cluster_ID *CT1_cluster_ID;
cluster_ID *CT2_cluster_ID;

void Read_CT1 (char *input) {
    FILE *fptr = fopen(input, "r");
    char Ori_Name[20];
    char Hm_Name[20];
    char col_name[20];
    int G_type;
    int Cond1_value[num_of_CT1_CMs];
    int Cond2_value[num_of_CT1_CMs];
    int total_columns, total_rows;

    printf("Reading CT1 table... ...");
    
    //reading head title
    total_columns = 2 * num_of_CT1_CMs + 3;
    CT1_cluster_ID = new cluster_ID[total_columns];
    
    for(int i=0; i< total_columns; i++) {
        fscanf(fptr, "%s", col_name);
        strcpy(CT1_cluster_ID[i].Col_Name, col_name);
    }
        
    num_of_CT1_Ligands = 0;
    num_of_CT1_Receptors = 0;
    while(fscanf(fptr,"%s", Ori_Name) != EOF) {
        fscanf(fptr,"%s", Hm_Name);
        fscanf(fptr,"\t%d", &G_type);
        for(int i=0; i<num_of_CT1_CMs; i++) {
            fscanf(fptr,"\t%d", &Cond1_value[i]);
        }
        for(int i=0; i<num_of_CT1_CMs; i++) {
            fscanf(fptr,"\t%d", &Cond2_value[i]);
        }
        if(G_type == 1) {
            num_of_CT1_Ligands++;
        }else if(G_type == 2) {
            num_of_CT1_Receptors++;
        }
    }

    total_rows = num_of_CT1_Ligands + num_of_CT1_Receptors;
    
    rewind(fptr);
    //table initialization
    
    CT1_Cond1_Ligand_exp_table = new Expression[num_of_CT1_Ligands];
    CT1_Cond2_Ligand_exp_table = new Expression[num_of_CT1_Ligands];
    CT1_Cond1_Receptor_exp_table = new Expression[num_of_CT1_Receptors];
    CT1_Cond2_Receptor_exp_table = new Expression[num_of_CT1_Receptors];
    
    for(int i=0; i<num_of_CT1_Ligands; i++) {
        CT1_Cond1_Ligand_exp_table[i].CM_exp = new int[num_of_CT1_CMs];
        CT1_Cond2_Ligand_exp_table[i].CM_exp = new int[num_of_CT1_CMs];
    }
    
    for(int i=0; i<num_of_CT1_Receptors; i++) {
        CT1_Cond1_Receptor_exp_table[i].CM_exp = new int[num_of_CT1_CMs];
        CT1_Cond2_Receptor_exp_table[i].CM_exp = new int[num_of_CT1_CMs];
    }
    
    for(int i=0; i<num_of_CT1_Ligands; i++) {
        for(int j=0; j<num_of_CT1_CMs; j++) {
            CT1_Cond1_Ligand_exp_table[i].CM_exp[j] = 0;
            CT1_Cond2_Ligand_exp_table[i].CM_exp[j] = 0;
        }
    }
    
    for(int i=0; i<num_of_CT1_Receptors; i++) {
        for(int j=0; j<num_of_CT1_CMs; j++) {
            CT1_Cond1_Receptor_exp_table[i].CM_exp[j] = 0;
            CT1_Cond2_Receptor_exp_table[i].CM_exp[j] = 0;
        }
    }
    
    
    //Read file and assign values to tables
    for(int i=0; i< total_columns; i++) {
        fscanf(fptr, "%s", col_name);
    }
    
    int index= 0;
    int index_receptor;
    while(fscanf(fptr,"%s", Ori_Name) != EOF) {
        if(index < num_of_CT1_Ligands) {
            strcpy(CT1_Cond1_Ligand_exp_table[index].original_name, Ori_Name);
            strcpy(CT1_Cond2_Ligand_exp_table[index].original_name, Ori_Name);
            fscanf(fptr,"%s", Hm_Name);
            strcpy(CT1_Cond1_Ligand_exp_table[index].gene_ID, Hm_Name);
            strcpy(CT1_Cond2_Ligand_exp_table[index].gene_ID, Hm_Name);
            fscanf(fptr,"\t%d", &G_type);
            for(int i=0; i<num_of_CT1_CMs; i++) {
                fscanf(fptr,"\t%d", &Cond1_value[i]);
                CT1_Cond1_Ligand_exp_table[index].CM_exp[i] = Cond1_value[i];
            }
            for(int i=0; i<num_of_CT1_CMs; i++) {
                fscanf(fptr,"\t%d", &Cond2_value[i]);
                CT1_Cond2_Ligand_exp_table[index].CM_exp[i] = Cond2_value[i];
            }
            index++;
        } else {
            index_receptor = index - num_of_CT1_Ligands;
            strcpy(CT1_Cond1_Receptor_exp_table[index_receptor].original_name, Ori_Name);
            strcpy(CT1_Cond2_Receptor_exp_table[index_receptor].original_name, Ori_Name);
            fscanf(fptr,"%s", Hm_Name);
            strcpy(CT1_Cond1_Receptor_exp_table[index_receptor].gene_ID, Hm_Name);
            strcpy(CT1_Cond2_Receptor_exp_table[index_receptor].gene_ID, Hm_Name);
            fscanf(fptr,"\t%d", &G_type);
            for(int i=0; i<num_of_CT1_CMs; i++) {
                fscanf(fptr,"\t%d", &Cond1_value[i]);
                CT1_Cond1_Receptor_exp_table[index_receptor].CM_exp[i] = Cond1_value[i];
            }
            for(int i=0; i<num_of_CT1_CMs; i++) {
                fscanf(fptr,"\t%d", &Cond2_value[i]);
                CT1_Cond2_Receptor_exp_table[index_receptor].CM_exp[i] = Cond2_value[i];
            }
            index++;
        }
    }
    
    fclose(fptr);
    
    printf(" ...done\n");
}

void Read_CT2 (char *input) {
    FILE *fptr = fopen(input, "r");
    char Ori_Name[20];
    char Hm_Name[20];
    char col_name[20];
    int G_type;
    int Cond1_value[num_of_CT2_CMs];
    int Cond2_value[num_of_CT2_CMs];
    int total_columns, total_rows;

    printf("Reading CT2 table... ...");
    
    //reading head title
    total_columns = 2 * num_of_CT2_CMs + 3;
    CT2_cluster_ID = new cluster_ID[total_columns];
    
    for(int i=0; i< total_columns; i++) {
        fscanf(fptr, "%s", col_name);
        strcpy(CT2_cluster_ID[i].Col_Name, col_name);
    }
        
    num_of_CT2_Ligands = 0;
    num_of_CT2_Receptors = 0;
    while(fscanf(fptr,"%s", Ori_Name) != EOF) {
        fscanf(fptr,"%s", Hm_Name);
        fscanf(fptr,"\t%d", &G_type);
        for(int i=0; i<num_of_CT2_CMs; i++) {
            fscanf(fptr,"\t%d", &Cond1_value[i]);
        }
        for(int i=0; i<num_of_CT2_CMs; i++) {
            fscanf(fptr,"\t%d", &Cond2_value[i]);
        }
        if(G_type == 1) {
            num_of_CT2_Ligands++;
        }else if(G_type == 2) {
            num_of_CT2_Receptors++;
        }
    }

    total_rows = num_of_CT2_Ligands + num_of_CT2_Receptors;
    
    rewind(fptr);
    //table initialization
    
    CT2_Cond1_Ligand_exp_table = new Expression[num_of_CT2_Ligands];
    CT2_Cond2_Ligand_exp_table = new Expression[num_of_CT2_Ligands];
    CT2_Cond1_Receptor_exp_table = new Expression[num_of_CT2_Receptors];
    CT2_Cond2_Receptor_exp_table = new Expression[num_of_CT2_Receptors];
    
    for(int i=0; i<num_of_CT2_Ligands; i++) {
        CT2_Cond1_Ligand_exp_table[i].CM_exp = new int[num_of_CT2_CMs];
        CT2_Cond2_Ligand_exp_table[i].CM_exp = new int[num_of_CT2_CMs];
    }
    
    for(int i=0; i<num_of_CT2_Receptors; i++) {
        CT2_Cond1_Receptor_exp_table[i].CM_exp = new int[num_of_CT2_CMs];
        CT2_Cond2_Receptor_exp_table[i].CM_exp = new int[num_of_CT2_CMs];
    }
    
    for(int i=0; i<num_of_CT2_Ligands; i++) {
        for(int j=0; j<num_of_CT2_CMs; j++) {
            CT2_Cond1_Ligand_exp_table[i].CM_exp[j] = 0;
            CT2_Cond2_Ligand_exp_table[i].CM_exp[j] = 0;
        }
    }
    
    for(int i=0; i<num_of_CT2_Receptors; i++) {
        for(int j=0; j<num_of_CT2_CMs; j++) {
            CT2_Cond1_Receptor_exp_table[i].CM_exp[j] = 0;
            CT2_Cond2_Receptor_exp_table[i].CM_exp[j] = 0;
        }
    }
    
    
    //Read file and assign values to tables
    for(int i=0; i< total_columns; i++) {
        fscanf(fptr, "%s", col_name);
    }
    
    int index= 0;
    int index_receptor;
    while(fscanf(fptr,"%s", Ori_Name) != EOF) {
        if(index < num_of_CT2_Ligands) {
            strcpy(CT2_Cond1_Ligand_exp_table[index].original_name, Ori_Name);
            strcpy(CT2_Cond2_Ligand_exp_table[index].original_name, Ori_Name);
            fscanf(fptr,"%s", Hm_Name);
            strcpy(CT2_Cond1_Ligand_exp_table[index].gene_ID, Hm_Name);
            strcpy(CT2_Cond2_Ligand_exp_table[index].gene_ID, Hm_Name);
            fscanf(fptr,"\t%d", &G_type);
            for(int i=0; i<num_of_CT2_CMs; i++) {
                fscanf(fptr,"\t%d", &Cond1_value[i]);
                CT2_Cond1_Ligand_exp_table[index].CM_exp[i] = Cond1_value[i];
            }
            for(int i=0; i<num_of_CT2_CMs; i++) {
                fscanf(fptr,"\t%d", &Cond2_value[i]);
                CT2_Cond2_Ligand_exp_table[index].CM_exp[i] = Cond2_value[i];
            }
            index++;
        } else {
            index_receptor = index - num_of_CT2_Ligands;
            strcpy(CT2_Cond1_Receptor_exp_table[index_receptor].original_name, Ori_Name);
            strcpy(CT2_Cond2_Receptor_exp_table[index_receptor].original_name, Ori_Name);
            fscanf(fptr,"%s", Hm_Name);
            strcpy(CT2_Cond1_Receptor_exp_table[index_receptor].gene_ID, Hm_Name);
            strcpy(CT2_Cond2_Receptor_exp_table[index_receptor].gene_ID, Hm_Name);
            fscanf(fptr,"\t%d", &G_type);
            for(int i=0; i<num_of_CT2_CMs; i++) {
                fscanf(fptr,"\t%d", &Cond1_value[i]);
                CT2_Cond1_Receptor_exp_table[index_receptor].CM_exp[i] = Cond1_value[i];
            }
            for(int i=0; i<num_of_CT2_CMs; i++) {
                fscanf(fptr,"\t%d", &Cond2_value[i]);
                CT2_Cond2_Receptor_exp_table[index_receptor].CM_exp[i] = Cond2_value[i];
            }
            index++;
        }
    }
    
    fclose(fptr);
    printf(" ...done\n");
}

void Read_LR_pair_DB (char *input) {
    FILE *fptr = fopen(input, "r");
    char LR_pair_name[40];
    
    printf("Reading known_L-R_pairs table... ...");
    
    num_of_LR_pairs = 0;
    while(fscanf(fptr,"%s", LR_pair_name) != EOF) {
        num_of_LR_pairs++;
    }
    
    rewind(fptr);
    //table initialization
    LR_table = new Pair[num_of_LR_pairs];
    
    int index = 0;
    while(fscanf(fptr,"%s", LR_pair_name) != EOF) {
        strcpy(LR_table[index].pair_name, LR_pair_name);
        strcpy(LR_table[index].pair_ori_name,"");
        LR_table[index].hit = 0;
        index++;
    }
    
    fclose(fptr);
    printf(" ...done\n");
}


int check_LR_pair_DB (char *pair_ID) {
    int hit_DB = 0;
    
    for(int i=0; i<num_of_LR_pairs; i++) {
        //if(strcmp(pair_ID, LR_table[i].pair_name) == 0 && LR_table[i].hit == 0) {
        if(strcmp(pair_ID, LR_table[i].pair_name) == 0) {
            //LR_table[i].hit = 1;
            hit_DB = 1;
        }
    }
    
    return hit_DB;
}


void pair_generation_1() {
    char lr_pair_for_test[60];
    char lr_pair_ori_name[60];
    int num_pair_cond1 = 0;
    int num_pair_cond2 = 0;
    int index_cond1 = 0;
    int index_cond2 = 0;
    
    FILE *fout1, *fout2, *fout3;
    
    fout1 = fopen("Cond1-specific_CT2-ligand-CT1-receptor_pairs.csv","w");
    fout2 = fopen("common_CT2-ligand-CT1-receptor_pairs.csv","w");
    fout3 = fopen("Cond2-specific_CT2-ligand-CT1-receptor_pairs.csv","w");
    
    printf("Generating CT2_ligand->CT1_receptor pairs... ...");
    
    for(int i=0; i<num_of_CT2_CMs; i++) {
        for(int k=0; k<num_of_CT1_CMs; k++) {
            
            for(int j=0; j<num_of_CT2_Ligands; j++) {
                if(CT2_Cond1_Ligand_exp_table[j].CM_exp[i] == 1) {
                    for(int l=0; l<num_of_CT1_Receptors; l++) {
                        if(CT1_Cond1_Receptor_exp_table[l].CM_exp[k] == 1) {
                            strcpy(lr_pair_for_test, "");
                            strcat(lr_pair_for_test, CT2_Cond1_Ligand_exp_table[j].gene_ID);
                            strcat(lr_pair_for_test, "_");
                            strcat(lr_pair_for_test, CT1_Cond1_Receptor_exp_table[l].gene_ID);
                            
                            if(check_LR_pair_DB(lr_pair_for_test) == 1) num_pair_cond1++;
                        }
                    }
                }
            }
            
            
            for(int j=0; j<num_of_CT2_Ligands; j++) {
                if(CT2_Cond2_Ligand_exp_table[j].CM_exp[i] == 1) {
                    for(int l=0; l<num_of_CT1_Receptors; l++) {
                        if(CT1_Cond2_Receptor_exp_table[l].CM_exp[k] == 1) {
                            strcpy(lr_pair_for_test, "");
                            strcat(lr_pair_for_test, CT2_Cond2_Ligand_exp_table[j].gene_ID);
                            strcat(lr_pair_for_test, "_");
                            strcat(lr_pair_for_test, CT1_Cond2_Receptor_exp_table[l].gene_ID);
                            
                            if(check_LR_pair_DB(lr_pair_for_test) == 1) num_pair_cond2++;
                        }
                    }
                }
            }
            
        }
    }
    
    Cond1_CT2_ligand_to_CT1_receptor_pair_list = new Pair[num_pair_cond1];
    Cond2_CT2_ligand_to_CT1_receptor_pair_list = new Pair[num_pair_cond2];
    
    for(int i=0; i<num_of_CT2_CMs; i++) {
        for(int k=0; k<num_of_CT1_CMs; k++) {
            for(int j=0; j<num_of_CT2_Ligands; j++) {
                if(CT2_Cond1_Ligand_exp_table[j].CM_exp[i] == 1) {
                    for(int l=0; l<num_of_CT1_Receptors; l++) {
                        if(CT1_Cond1_Receptor_exp_table[l].CM_exp[k] == 1) {
                            strcpy(lr_pair_for_test, "");
                            strcat(lr_pair_for_test, CT2_Cond1_Ligand_exp_table[j].gene_ID);
                            strcat(lr_pair_for_test, "_");
                            strcat(lr_pair_for_test, CT1_Cond1_Receptor_exp_table[l].gene_ID);
                            
                            strcpy(lr_pair_ori_name, "");
                            strcat(lr_pair_ori_name, CT2_Cond1_Ligand_exp_table[j].original_name);
                            strcat(lr_pair_ori_name, "_");
                            strcat(lr_pair_ori_name, CT1_Cond1_Receptor_exp_table[l].original_name);
                            
                            if(check_LR_pair_DB(lr_pair_for_test) == 1) {
                                
                                strcpy(Cond1_CT2_ligand_to_CT1_receptor_pair_list[index_cond1].pair_name,lr_pair_for_test);
                                strcpy(Cond1_CT2_ligand_to_CT1_receptor_pair_list[index_cond1].pair_ori_name,lr_pair_ori_name);
                                strcpy(Cond1_CT2_ligand_to_CT1_receptor_pair_list[index_cond1].cluster_ligand, CT2_cluster_ID[i+3].Col_Name);
                                strcpy(Cond1_CT2_ligand_to_CT1_receptor_pair_list[index_cond1].cluster_receptor, CT1_cluster_ID[k+3].Col_Name);
                                Cond1_CT2_ligand_to_CT1_receptor_pair_list[index_cond1].cluster_ID = i * num_of_CT1_CMs + k;
                                index_cond1++;
                            }
                        }
                    }
                }
            }

            for(int j=0; j<num_of_CT2_Ligands; j++) {
                if(CT2_Cond2_Ligand_exp_table[j].CM_exp[i] == 1) {
                    for(int l=0; l<num_of_CT1_Receptors; l++) {
                        if(CT1_Cond2_Receptor_exp_table[l].CM_exp[k] == 1) {
                            strcpy(lr_pair_for_test, "");
                            strcat(lr_pair_for_test, CT2_Cond2_Ligand_exp_table[j].gene_ID);
                            strcat(lr_pair_for_test, "_");
                            strcat(lr_pair_for_test, CT1_Cond2_Receptor_exp_table[l].gene_ID);
                            
                            strcpy(lr_pair_ori_name, "");
                            strcat(lr_pair_ori_name, CT2_Cond2_Ligand_exp_table[j].original_name);
                            strcat(lr_pair_ori_name, "_");
                            strcat(lr_pair_ori_name, CT1_Cond2_Receptor_exp_table[l].original_name);
                            
                            if(check_LR_pair_DB(lr_pair_for_test) == 1) {
                                
                                strcpy(Cond2_CT2_ligand_to_CT1_receptor_pair_list[index_cond2].pair_name,lr_pair_for_test);
                                strcpy(Cond2_CT2_ligand_to_CT1_receptor_pair_list[index_cond2].pair_ori_name,lr_pair_ori_name);
                                strcpy(Cond2_CT2_ligand_to_CT1_receptor_pair_list[index_cond2].cluster_ligand, CT2_cluster_ID[i+3].Col_Name);
                                strcpy(Cond2_CT2_ligand_to_CT1_receptor_pair_list[index_cond2].cluster_receptor, CT1_cluster_ID[k+3].Col_Name);
                                Cond2_CT2_ligand_to_CT1_receptor_pair_list[index_cond2].cluster_ID = i * num_of_CT1_CMs + k;
                                index_cond2++;
                            }
                        }
                    }
                }
            }
            
        }
    }
    printf(" ...done\n");
    
    printf("Identifying condition-specific and common pairs... ...");
    // Compare two L-R_pair_lists between condition 1 and 2
    for(int i=0; i<num_pair_cond1; i++) Cond1_CT2_ligand_to_CT1_receptor_pair_list[i].hit = 0;
    for(int i=0; i<num_pair_cond2; i++) Cond2_CT2_ligand_to_CT1_receptor_pair_list[i].hit = 0;
    
    for(int i=0; i<num_pair_cond1; i++) {
        for(int j=0; j<num_pair_cond2; j++) {
            if(Cond1_CT2_ligand_to_CT1_receptor_pair_list[i].cluster_ID == Cond2_CT2_ligand_to_CT1_receptor_pair_list[j].cluster_ID) {
                if(strcmp(Cond1_CT2_ligand_to_CT1_receptor_pair_list[i].pair_ori_name, Cond2_CT2_ligand_to_CT1_receptor_pair_list[j].pair_ori_name) == 0) {
                    Cond1_CT2_ligand_to_CT1_receptor_pair_list[i].hit = 1;
                    Cond2_CT2_ligand_to_CT1_receptor_pair_list[j].hit = 1;
                }
            }
        }
    }
    
    fprintf(fout1, "Cluster pair ID,CT2,pair type,CT1,human L-R pair,L-R pair\n");
    fprintf(fout2, "Cluster pair ID,CT2,pair type,CT1,human L-R pair,L-R pair\n");
    fprintf(fout3, "Cluster pair ID,CT2,pair type,CT1,human L-R pair,L-R pair\n");
    
    // Calculate and print condition-specific and common pairs
    for (int i=0; i<num_pair_cond1; i++) {
        if(Cond1_CT2_ligand_to_CT1_receptor_pair_list[i].hit == 0) {
            fprintf(fout1, "%d", Cond1_CT2_ligand_to_CT1_receptor_pair_list[i].cluster_ID);
            fprintf(fout1, ",%s", Cond1_CT2_ligand_to_CT1_receptor_pair_list[i].cluster_ligand);
            fprintf(fout1, ",L->R");
            fprintf(fout1, ",%s", Cond1_CT2_ligand_to_CT1_receptor_pair_list[i].cluster_receptor);
            fprintf(fout1, ",%s", Cond1_CT2_ligand_to_CT1_receptor_pair_list[i].pair_name);
            fprintf(fout1, ",%s\n", Cond1_CT2_ligand_to_CT1_receptor_pair_list[i].pair_ori_name);
        } else {
            fprintf(fout2, "%d", Cond1_CT2_ligand_to_CT1_receptor_pair_list[i].cluster_ID);
            fprintf(fout2, ",%s", Cond1_CT2_ligand_to_CT1_receptor_pair_list[i].cluster_ligand);
            fprintf(fout2, ",L->R");
            fprintf(fout2, ",%s", Cond1_CT2_ligand_to_CT1_receptor_pair_list[i].cluster_receptor);
            fprintf(fout2, ",%s", Cond1_CT2_ligand_to_CT1_receptor_pair_list[i].pair_name);
            fprintf(fout2, ",%s\n", Cond1_CT2_ligand_to_CT1_receptor_pair_list[i].pair_ori_name);
        }
    }
    
    for (int i=0; i<num_pair_cond2; i++) {
        if(Cond2_CT2_ligand_to_CT1_receptor_pair_list[i].hit == 0) {
            fprintf(fout3, "%d", Cond2_CT2_ligand_to_CT1_receptor_pair_list[i].cluster_ID);
            fprintf(fout3, ",%s", Cond2_CT2_ligand_to_CT1_receptor_pair_list[i].cluster_ligand);
            fprintf(fout3, ",L->R");
            fprintf(fout3, ",%s", Cond2_CT2_ligand_to_CT1_receptor_pair_list[i].cluster_receptor);
            fprintf(fout3, ",%s", Cond2_CT2_ligand_to_CT1_receptor_pair_list[i].pair_name);
            fprintf(fout3, ",%s\n", Cond2_CT2_ligand_to_CT1_receptor_pair_list[i].pair_ori_name);
        }
    }
    
    printf(" ...done\n");
    
    fclose(fout1);
    fclose(fout2);
    fclose(fout3);
}

void pair_generation_2() {
    char lr_pair_for_test[60];
    char lr_pair_ori_name[60];
    int num_pair_cond1 = 0;
    int num_pair_cond2 = 0;
    int index_cond1 = 0;
    int index_cond2 = 0;
    
    FILE *fout1, *fout2, *fout3;
    
    fout1 = fopen("Cond1-specific_CT1-ligand-CT2-receptor_pairs.csv","w");
    fout2 = fopen("common_CT1-ligand-CT2-receptor_pairs.csv","w");
    fout3 = fopen("Cond2-specific_CT1-ligand-CT2-receptor_pairs.csv","w");
    
    printf("Generating CT1_ligand->CT2_receptor pairs... ...");
    
    for(int i=0; i<num_of_CT1_CMs; i++) {
        for(int k=0; k<num_of_CT2_CMs; k++) {
            
            for(int j=0; j<num_of_CT1_Ligands; j++) {
                if(CT1_Cond1_Ligand_exp_table[j].CM_exp[i] == 1) {
                    for(int l=0; l<num_of_CT2_Receptors; l++) {
                        if(CT2_Cond1_Receptor_exp_table[l].CM_exp[k] == 1) {
                            strcpy(lr_pair_for_test, "");
                            strcat(lr_pair_for_test, CT1_Cond1_Ligand_exp_table[j].gene_ID);
                            strcat(lr_pair_for_test, "_");
                            strcat(lr_pair_for_test, CT2_Cond1_Receptor_exp_table[l].gene_ID);
                            
                            if(check_LR_pair_DB(lr_pair_for_test) == 1) num_pair_cond1++;
                        }
                    }
                }
            }
            
            
            for(int j=0; j<num_of_CT1_Ligands; j++) {
                if(CT1_Cond2_Ligand_exp_table[j].CM_exp[i] == 1) {
                    for(int l=0; l<num_of_CT2_Receptors; l++) {
                        if(CT2_Cond2_Receptor_exp_table[l].CM_exp[k] == 1) {
                            strcpy(lr_pair_for_test, "");
                            strcat(lr_pair_for_test, CT1_Cond2_Ligand_exp_table[j].gene_ID);
                            strcat(lr_pair_for_test, "_");
                            strcat(lr_pair_for_test, CT2_Cond2_Receptor_exp_table[l].gene_ID);
                            
                            if(check_LR_pair_DB(lr_pair_for_test) == 1) num_pair_cond2++;
                        }
                    }
                }
            }
            
        }
    }
    
    Cond1_CT1_ligand_to_CT2_receptor_pair_list = new Pair[num_pair_cond1];
    Cond2_CT1_ligand_to_CT2_receptor_pair_list = new Pair[num_pair_cond2];
    
    for(int i=0; i<num_of_CT1_CMs; i++) {
        for(int k=0; k<num_of_CT2_CMs; k++) {
            for(int j=0; j<num_of_CT1_Ligands; j++) {
                if(CT1_Cond1_Ligand_exp_table[j].CM_exp[i] == 1) {
                    for(int l=0; l<num_of_CT2_Receptors; l++) {
                        if(CT2_Cond1_Receptor_exp_table[l].CM_exp[k] == 1) {
                            strcpy(lr_pair_for_test, "");
                            strcat(lr_pair_for_test, CT1_Cond1_Ligand_exp_table[j].gene_ID);
                            strcat(lr_pair_for_test, "_");
                            strcat(lr_pair_for_test, CT2_Cond1_Receptor_exp_table[l].gene_ID);
                            
                            strcpy(lr_pair_ori_name, "");
                            strcat(lr_pair_ori_name, CT1_Cond1_Ligand_exp_table[j].original_name);
                            strcat(lr_pair_ori_name, "_");
                            strcat(lr_pair_ori_name, CT2_Cond1_Receptor_exp_table[l].original_name);
                            
                            if(check_LR_pair_DB(lr_pair_for_test) == 1) {
                                
                                strcpy(Cond1_CT1_ligand_to_CT2_receptor_pair_list[index_cond1].pair_name,lr_pair_for_test);
                                strcpy(Cond1_CT1_ligand_to_CT2_receptor_pair_list[index_cond1].pair_ori_name,lr_pair_ori_name);
                                strcpy(Cond1_CT1_ligand_to_CT2_receptor_pair_list[index_cond1].cluster_ligand, CT1_cluster_ID[i+3].Col_Name);
                                strcpy(Cond1_CT1_ligand_to_CT2_receptor_pair_list[index_cond1].cluster_receptor, CT2_cluster_ID[k+3].Col_Name);
                                Cond1_CT1_ligand_to_CT2_receptor_pair_list[index_cond1].cluster_ID = i * num_of_CT2_CMs + k;
                                index_cond1++;
                            }
                        }
                    }
                }
            }

            for(int j=0; j<num_of_CT1_Ligands; j++) {
                if(CT1_Cond2_Ligand_exp_table[j].CM_exp[i] == 1) {
                    for(int l=0; l<num_of_CT2_Receptors; l++) {
                        if(CT2_Cond2_Receptor_exp_table[l].CM_exp[k] == 1) {
                            strcpy(lr_pair_for_test, "");
                            strcat(lr_pair_for_test, CT1_Cond2_Ligand_exp_table[j].gene_ID);
                            strcat(lr_pair_for_test, "_");
                            strcat(lr_pair_for_test, CT2_Cond2_Receptor_exp_table[l].gene_ID);
                            
                            strcpy(lr_pair_ori_name, "");
                            strcat(lr_pair_ori_name, CT1_Cond2_Ligand_exp_table[j].original_name);
                            strcat(lr_pair_ori_name, "_");
                            strcat(lr_pair_ori_name, CT2_Cond2_Receptor_exp_table[l].original_name);
                            
                            if(check_LR_pair_DB(lr_pair_for_test) == 1) {
                                
                                strcpy(Cond2_CT1_ligand_to_CT2_receptor_pair_list[index_cond2].pair_name,lr_pair_for_test);
                                strcpy(Cond2_CT1_ligand_to_CT2_receptor_pair_list[index_cond2].pair_ori_name,lr_pair_ori_name);
                                strcpy(Cond2_CT1_ligand_to_CT2_receptor_pair_list[index_cond2].cluster_ligand, CT1_cluster_ID[i+3].Col_Name);
                                strcpy(Cond2_CT1_ligand_to_CT2_receptor_pair_list[index_cond2].cluster_receptor, CT2_cluster_ID[k+3].Col_Name);
                                Cond2_CT1_ligand_to_CT2_receptor_pair_list[index_cond2].cluster_ID = i * num_of_CT2_CMs + k;
                                index_cond2++;
                            }
                        }
                    }
                }
            }
            
        }
    }
    printf(" ...done\n");
    
    printf("Identifying condition-specific and common pairs... ...");
    // Compare two L-R_pair_lists between condition 1 and 2
    for(int i=0; i<num_pair_cond1; i++) Cond1_CT1_ligand_to_CT2_receptor_pair_list[i].hit = 0;
    for(int i=0; i<num_pair_cond2; i++) Cond2_CT1_ligand_to_CT2_receptor_pair_list[i].hit = 0;
    
    for(int i=0; i<num_pair_cond1; i++) {
        for(int j=0; j<num_pair_cond2; j++) {
            if(Cond1_CT1_ligand_to_CT2_receptor_pair_list[i].cluster_ID == Cond2_CT1_ligand_to_CT2_receptor_pair_list[j].cluster_ID) {
                if(strcmp(Cond1_CT1_ligand_to_CT2_receptor_pair_list[i].pair_ori_name, Cond2_CT1_ligand_to_CT2_receptor_pair_list[j].pair_ori_name) == 0) {
                    Cond1_CT1_ligand_to_CT2_receptor_pair_list[i].hit = 1;
                    Cond2_CT1_ligand_to_CT2_receptor_pair_list[j].hit = 1;
                }
            }
        }
    }
    
    fprintf(fout1, "Cluster pair ID,CT1,pair type,CT2,human L-R pair,L-R pair\n");
    fprintf(fout2, "Cluster pair ID,CT1,pair type,CT2,human L-R pair,L-R pair\n");
    fprintf(fout3, "Cluster pair ID,CT1,pair type,CT2,human L-R pair,L-R pair\n");
    
    // Calculate and print condition-specific and common pairs
    for (int i=0; i<num_pair_cond1; i++) { // Cond1-specific set
        if(Cond1_CT1_ligand_to_CT2_receptor_pair_list[i].hit == 0) {
            fprintf(fout1, "%d", Cond1_CT1_ligand_to_CT2_receptor_pair_list[i].cluster_ID);
            fprintf(fout1, ",%s", Cond1_CT1_ligand_to_CT2_receptor_pair_list[i].cluster_ligand);
            fprintf(fout1, ",L->R");
            fprintf(fout1, ",%s", Cond1_CT1_ligand_to_CT2_receptor_pair_list[i].cluster_receptor);
            fprintf(fout1, ",%s", Cond1_CT1_ligand_to_CT2_receptor_pair_list[i].pair_name);
            fprintf(fout1, ",%s\n", Cond1_CT1_ligand_to_CT2_receptor_pair_list[i].pair_ori_name);
        } else { // Common set
            fprintf(fout2, "%d", Cond1_CT1_ligand_to_CT2_receptor_pair_list[i].cluster_ID);
            fprintf(fout2, ",%s", Cond1_CT1_ligand_to_CT2_receptor_pair_list[i].cluster_ligand);
            fprintf(fout2, ",L->R");
            fprintf(fout2, ",%s", Cond1_CT1_ligand_to_CT2_receptor_pair_list[i].cluster_receptor);
            fprintf(fout2, ",%s", Cond1_CT1_ligand_to_CT2_receptor_pair_list[i].pair_name);
            fprintf(fout2, ",%s\n", Cond1_CT1_ligand_to_CT2_receptor_pair_list[i].pair_ori_name);
        }
    }
    
    for (int i=0; i<num_pair_cond2; i++) { //Cond2-speicific set
        if(Cond2_CT1_ligand_to_CT2_receptor_pair_list[i].hit == 0) {
            fprintf(fout3, "%d", Cond2_CT1_ligand_to_CT2_receptor_pair_list[i].cluster_ID);
            fprintf(fout3, ",%s", Cond2_CT1_ligand_to_CT2_receptor_pair_list[i].cluster_ligand);
            fprintf(fout3, ",L->R");
            fprintf(fout3, ",%s", Cond2_CT1_ligand_to_CT2_receptor_pair_list[i].cluster_receptor);
            fprintf(fout3, ",%s", Cond2_CT1_ligand_to_CT2_receptor_pair_list[i].pair_name);
            fprintf(fout3, ",%s\n", Cond2_CT1_ligand_to_CT2_receptor_pair_list[i].pair_ori_name);
        }
    }
    
    printf(" ...done\n");
    
    fclose(fout1);
    fclose(fout2);
    fclose(fout3);
}

void debug() {
    FILE *fout1;
    
    fout1 = fopen("debug_out.txt","w");
    
    for(int i=0; i<num_of_CT2_Receptors; i++) {
        fprintf(fout1, "%s, %s",CT2_Cond2_Receptor_exp_table[i].original_name, CT2_Cond2_Receptor_exp_table[i].gene_ID);
        for(int j=0; j<num_of_CT2_CMs; j++) {
            fprintf(fout1, ", %d", CT2_Cond2_Receptor_exp_table[i].CM_exp[j]);
        }
        fprintf(fout1, "\n");
    }
    fclose(fout1);
}

int main(int argc, char* argv[]) {
    char input_file1[100];   //CT1 gene list
    char input_file2[100];   //CT2 gene list
    char input_file3[100];   //ligand-receptor pair from database
    
    
    if (argc != 6) {
        printf("\nUsage: LR_generator No_of_CT1_clusters No_of_CT2_clusters CT1_table_file CT2_table_file L-R_pair_table\n");
    } else {
        
        num_of_CT1_CMs = atoi(argv[1]);
        num_of_CT2_CMs = atoi(argv[2]);
        strcpy(input_file1, argv[3]);
        strcpy(input_file2, argv[4]);
        strcpy(input_file3, argv[5]);
    
        FILE *fptr1 = fopen(input_file1, "r");
        FILE *fptr2 = fopen(input_file2, "r");
        FILE *fptr3 = fopen(input_file3, "r");
        
        if(fptr1 == NULL || fptr2 == NULL || fptr3 == NULL) {
            
            printf("\nCan't find the input file. Please check the inupt file again!\n\n");
            
        } else {
    
            fclose (fptr1);
            fclose (fptr2);
            fclose (fptr3);
    
            Read_CT1(input_file1);
            Read_CT2(input_file2);
            Read_LR_pair_DB(input_file3);
            //debug();
    
            pair_generation_1();
            pair_generation_2();
        }
    }
    
    return 0;
}
