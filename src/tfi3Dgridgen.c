#include "tfi2Dgridgen.h"
#include "tfi3Dgridgen.h"


void generate_EMC3_neutral_3Dgrid_TFI(ThreeDimGrid* grid3d_neu,
                                      const ThreeDimGrid* grid3d_plasma,
                                      DLListNode* inner_tgt_head, DLListNode* outer_tgt_head, DLListNode* neu_top_head,
                                      const int n_neu_left_distrb, const double* neu_left_distrb,
                                      ode_function* func, ode_solver* solver)
{

//1. Check the inputs
  if(!grid3d_neu || !grid3d_plasma)
  {
    fprintf(stderr, "Error: Empty 3D grid in generate_EMC3_neutral_3Dgrid_TFI.\n");
    exit(EXIT_FAILURE);
  }

  if(!inner_tgt_head || !outer_tgt_head || !neu_top_head)
  {
    fprintf(stderr, "Error: Empty boundary in generate_EMC3_neutral_3Dgrid_TFI.\n");
    exit(EXIT_FAILURE);
  }

  if(n_neu_left_distrb<2 || !neu_left_distrb)
  {
    fprintf(stderr, "Error: Unexpected inputs of neutral grid distribution in generate_EMC3_neutral_3Dgrid_TFI.\n");
    exit(EXIT_FAILURE);
  }

  if(grid3d_neu->npol != grid3d_plasma->npol ||
     grid3d_neu->ntor != grid3d_plasma->ntor)
  {
    fprintf(stderr, "Error: Poloidal or Toroidal size of 3D grid between neutral and plasma are not identical.\n");
    exit(EXIT_FAILURE);
  }

  if(n_neu_left_distrb != grid3d_neu->nrad)
  {
    fprintf(stderr, "Error: The size of radial distribution is not identical with neutral 3D grid.\n");
    exit(EXIT_FAILURE);
  }

  bool is_SOL=false;
  bool is_PFR=false;

  DLListNode *neu_top_tail = get_DLList_tailnode(neu_top_head);
  DLListNode *inner_tgt_tail = get_DLList_tailnode(inner_tgt_head);
  DLListNode *outer_tgt_tail = get_DLList_tailnode(outer_tgt_head);

  if(fabs(neu_top_head->r - inner_tgt_head->r) < EPSILON_12 && 
     fabs(neu_top_head->z - inner_tgt_head->z) < EPSILON_12 &&
     fabs(neu_top_tail->r - outer_tgt_head->r) < EPSILON_12 &&
     fabs(neu_top_tail->z - outer_tgt_head->z) < EPSILON_12 )
  {
  #ifdef DEBUG
    printf("This is neutral grid for PFR region.\n");
  #endif
    is_PFR = true;
  }
  else if (fabs(neu_top_head->r - inner_tgt_tail->r) < EPSILON_12 && 
           fabs(neu_top_head->z - inner_tgt_tail->z) < EPSILON_12 &&
           fabs(neu_top_tail->r - outer_tgt_tail->r) < EPSILON_12 &&
           fabs(neu_top_tail->z - outer_tgt_tail->z) < EPSILON_12 )
  {
  #ifdef DEBUG
    printf("This is neutral grid for SOL region.\n");
  #endif
    is_SOL = true;
  }
  else
  {
    fprintf(stderr, "Error: The boundary of neutral grid is not connected with inner target or outer target in generate_EMC3_neutral_3Dgrid_TFI.\n");
    exit(EXIT_FAILURE);
  }
  
  //Extract phi values from plasma 3D grid
  int nphi=grid3d_plasma->ntor;

  double* phi=malloc(nphi*sizeof(double));
  for(int i=0; i<nphi; i++)
  {
    phi[i]=get_phi_3Dgrid(grid3d_plasma, 0,0,i);
  }


//2. Generate the neutral 3D grid

  //Temperoray 2d grid to store neu 2d grid,
  TwoDimGrid* grid2d_tmp=create_2Dgrid_poloidal_major(grid3d_plasma->npol-(nphi-1),n_neu_left_distrb);
  int npol = grid3d_plasma->npol;
  int nrad = grid3d_plasma->nrad;

  for(int k=0;k<nphi;k++) //in each toridal slice which is the RZ plance at fixed phi
  {
  // 2.1 Build the necessary inputs for TFI part
    if(k!=0) reset_2Dgrid(grid2d_tmp);
  //Extract the DDList for bottom boundary of neutral grid
    double r_tmp = get_r_3Dgrid(grid3d_plasma, 0, nrad-1, k);
    double z_tmp = get_z_3Dgrid(grid3d_plasma, 0, nrad-1, k);
    double phi_tmp = get_phi_3Dgrid(grid3d_plasma, 0, nrad-1, k);

    DLListNode* neu_bottom_head = create_DLListNode(r_tmp, z_tmp);
    DLListNode* tail_tmp = neu_bottom_head;
    for(int i=1;i<npol;i++)
    {
      add_DLListnode_at_tail(&tail_tmp, 
                              get_r_3Dgrid(grid3d_plasma, i, nrad-1, k), 
                              get_z_3Dgrid(grid3d_plasma, i, nrad-1, k));
    }

  //neu_bottom should NOT intersect with neu_top (neu_top_head)
    if(has_intersection_DLList(neu_top_head, neu_bottom_head)==0)
    {
      fprintf(stderr, "Error: In the toridal %d slice (phi=%lf), the bottom boundary intersect with top boundary.\n",k,phi[k]);
      #ifdef DEBUG
      char name[32];
      sprintf(name,"DEBUG_neu_bottom_%d",k);
      write_DLList(neu_bottom_head,name);
      sprintf(name,"DEBUG_neu_top_%d",k);
      write_DLList(neu_top_head,name);
      fprintf(stderr, "Please check the file DEBUG_neu_bottom_%d and DEBUG_neu_top_%d.\n",k,k);
      #endif
      exit(EXIT_FAILURE);
    }

  //Build the necessary boundarirs and cut botoom,left and right
    double left_intsct_r;
    double left_intsct_z;
    double right_intsct_r;
    double right_intsct_z;

    DLListNode* neu_left_head = copy_DLList(inner_tgt_head);
    DLListNode* neu_right_head = copy_DLList(outer_tgt_head);

    if(insert_intersections_DLList(neu_bottom_head,neu_left_head, &left_intsct_r, &left_intsct_z)!=0)
    {
      fprintf(stderr, "Unexpected Error: The bottom bnd do not intersect with inner target.\n");
      exit(EXIT_FAILURE);
    }

    if(insert_intersections_DLList(neu_bottom_head,neu_right_head, &right_intsct_r, &right_intsct_z)!=0)
    {
      fprintf(stderr, "Unexpected Error: The bottom bnd do not intersect with outer target.\n");
      exit(EXIT_FAILURE);
    }
    //Check the intersections
    //TODO

    // For SOL neutral, cut left and right
    if(is_SOL && (!is_PFR))
    {
      cut_DLList_before_point(&neu_left_head,left_intsct_r,left_intsct_z);
      cut_DLList_before_point(&neu_right_head,right_intsct_r,right_intsct_z);
    }
    // For PFR neutral, cut left and right
    else if(!is_SOL && is_PFR)
    {
      cut_DLList_after_point(neu_left_head,left_intsct_r,left_intsct_z);
      cut_DLList_after_point(neu_right_head,right_intsct_r,right_intsct_z);
      //reverse the direction
      reverse_DLList(&neu_left_head);
      reverse_DLList(&neu_right_head);
    }
    else
    {
      fprintf(stderr, "Unexpected Error: Unrecognize the region for neutral grid.\n");
      exit(EXIT_FAILURE);
    }
    // Cut bottom
    cut_DLList_before_point(&neu_bottom_head, left_intsct_r,left_intsct_z);
    cut_DLList_after_point(neu_bottom_head, right_intsct_r,right_intsct_z);

    #ifdef DEBUG
    char name[32];
    sprintf(name,"NEU_GIRD2D_BND_BOTTOM%d",k);
    write_DLList(neu_bottom_head,name);

    sprintf(name,"NEU_GIRD2D_BND_LEFT%d",k);
    write_DLList(neu_left_head,name);

    sprintf(name,"NEU_GIRD2D_BND_RIGHT%d",k);
    write_DLList(neu_right_head,name);

    #endif

    Curve* neu_top_curve=convert_ddl_to_curve(neu_top_head);
    Curve* neu_bottom_curve=convert_ddl_to_curve(neu_bottom_head);
    Curve* neu_left_curve=convert_ddl_to_curve(neu_left_head);
    Curve* neu_right_curve=convert_ddl_to_curve(neu_right_head);

    //check again
    if(neu_bottom_curve->n_point!=grid2d_tmp->npol)
    {
      fprintf(stderr, "Unexpected Error: The poloidal size of grid2d tmp is not same size of bottom boundary.\n");
      exit(EXIT_FAILURE);
    }

    //bottom distribution is same with bottom curve
    int nbottom=neu_bottom_curve->n_point;
    double* distrb_b=malloc((nbottom)*sizeof(double));
    double len_tmp = total_length_curve(neu_bottom_curve);

    for(int i=0; i<nbottom; i++)
    {
      distrb_b[i]=length_curve(neu_bottom_curve, i+1)/len_tmp;
    }
    
    distrb_b[0]=0.0;
    distrb_b[nbottom-1]=1.0;

  // 2.2 Generate the 2D neutral grid by TFI
    generate_2Dgrid_default_TFI(grid2d_tmp,
                                neu_bottom_curve, neu_top_curve, neu_left_curve, neu_right_curve,
                                distrb_b, nbottom,
                                distrb_b, nbottom,
                                neu_left_distrb,n_neu_left_distrb, 
                                neu_left_distrb,n_neu_left_distrb);
                                
    optimized_neu_2Dgrid(grid2d_tmp);

    #ifdef DEBUG
    sprintf(name,"NEU_GIRD2D_BASE%d",k);
    write_2Dgrid(grid2d_tmp,name);
    #endif


  // 2.3 Extend at the inner and outer target along magnetic field line.
  /*
  *  This part is very tricky. This is neutral grid, it is not mandatory the neutral grid along the magnetic field line.
  *  However, in order the ensure the neutral grid to match the size of 3D plasma grid we use similar way as plasme grid,
  *  The grid out of inner and outer target is determined by magnetic filed line.
  *  It can be artifical, but this method can avoid intersections because magnetic field line never intersect with each other.
  */

    expand_target_EMC3_2Dgrid_default(grid2d_tmp, func, solver, phi_tmp, nphi, phi);
    #ifdef DEBUG
    sprintf(name,"NEU_GIRD2D_EXPTGT%d",k);
    write_2Dgrid(grid2d_tmp,name);
    #endif

  // 2.4 assign the 3D grid based on the 2D grid
    assign_2D_to_3D_tor_slice(grid2d_tmp, grid3d_neu, k, phi_tmp);

  // 2.5 free uncessary paramters
    free_DLList(neu_bottom_head);
    free_DLList(neu_left_head);
    free_DLList(neu_right_head);
    free_curve(neu_top_curve);
    free_curve(neu_bottom_curve);
    free_curve(neu_left_curve);
    free_curve(neu_right_curve);
    free(distrb_b);
  }
//free parameters
  free(phi);
  free_2Dgrid(grid2d_tmp);
}

