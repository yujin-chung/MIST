/* MIST copyright 2016 by Yujin Chung and Jody Hey */


#include <iostream>
#include <math.h>
#include "histograms.hpp"








void hist::initiateX(unsigned int gridsize, double priorMax, double priorMin)
{
  if(priorMax == 0)
    {
      x.resize(1,1);
      x(0,0) = 0.0;
    }
  else
    {
      x.resize(1,gridsize);
      for(unsigned int i=0; i<gridsize; i++)
	{
	  x(0,i) =  priorMin + ((i + 0.5) * (priorMax - priorMin)) / gridsize;
	}
    }
  return;
}



unsigned int hist::find_binID(double val, double numBins, double priorMax, double priorMin)
{
  unsigned int id;
  // find the id of bin where "val" falls
  id = static_cast<unsigned int> ((val-priorMin)/((priorMax - priorMin) / numBins));


  return id;
}




void marginal::initiate_Points2Integrate_largeBinWidth()
{
	for(unsigned int id =0; id<nPara; id++)
	{
		unsigned int binWidth = gridsize/nJointPoints;
		binWidth_large.at(id) = binWidth * binWidth_small.at(id); // binWidth_large = 0 if priorMax = 0
		if(priorsMax.at(id) == 0.0)
		  {
		    points2Integrate.at(id).push_back(0.0);
		  }
		else
		  {
		    for(unsigned int i=0; i<nJointPoints; i++)
		      points2Integrate.at(id).push_back(hp.at(id).getX(binWidth/2 + i*binWidth));
		  }
	}
	return;
}

void marginal::initiate(IM im)
{
  // gridsize = 10;//6;//50;//100;
  // nJointPoints =2;//2;// 3;//4;
  nPara = 6; // 3 population sizes, 2 migration rates, 1 splitting time.
  
  for(unsigned int i=0; i<nPara; i++)
    {
      if(i < 3)
	priorsMax.push_back(im.get_popSizeMax());
      else if(i==nPara-1)
	priorsMax.push_back(im.get_splittingTimeMax());
      else
	priorsMax.push_back(im.get_migRateMax());
    }
  hp.resize(nPara);
  // jointDensity.resize(0,nPara+1);
  points2Integrate.resize(nPara);
  binWidth_small.resize(nPara);
  binWidth_large.resize(nPara);
  
  unsigned int zeroPara = 0;
  for(unsigned int i=0; i<nPara;i++)
    if(priorsMax.at(i) == 0.0)
      zeroPara++;
  nZeroPara = zeroPara;

  unsigned int nrow = pow(nJointPoints,nPara-zeroPara);
  jointDensity.resize(nrow,nPara+1);
  jointDensity.setZero();
  
  return;
}

void marginal::initiate_forDE(IM im, unsigned int numPara, std::vector<double> paraMax)
{
  nPara = numPara;
  
  for(unsigned int i=0; i<nPara; i++)
    {
      priorsMax.push_back(paraMax.at(i));
    }
  hp.resize(nPara);
  //points2Integrate.resize(nPara);
  binWidth_small.resize(nPara);
  // binWidth_large.resize(nPara);

  gridsize = 1000;
  initiate_hist_forDE();

  /* 
  unsigned int zeroPara = 0;
  for(unsigned int i=0; i<nPara;i++)
    if(priorsMax.at(i) == 0.0)
      zeroPara++;
  nZeroPara = zeroPara;
  */

  /*
  unsigned int nrow = pow(nJointPoints,nPara-zeroPara);
  jointDensity.resize(nrow,nPara+1);
  jointDensity.setZero();
  */
  
  return;
}

void marginal::initiate_hist_forDE()
{
  double priorMin = 0.0;
  for(unsigned int i=0; i<nPara; i++)
    {
      hp.at(i).initiateX(gridsize, priorsMax.at(i),priorMin);
      hp.at(i).initiateY(gridsize);
      hp.at(i).initiate_freq(gridsize);
      
      double width = (priorsMax.at(i)-priorMin)/gridsize;
      binWidth_small.at(i) = width; // width =0 if priorMax = 0
    }
  return;  
}

void marginal::initiate_histXvalues_binWidth_small()
{
	double priorMin = 0.0;
	for(unsigned int i=0; i<nPara; i++)
	{
	  hp.at(i).initiateX(gridsize, priorsMax.at(i),priorMin);
	  
	  double width = (priorsMax.at(i)-priorMin)/gridsize;
	  binWidth_small.at(i) = width; // width =0 if priorMax = 0
	}
	return;
}

void marginal::add(Eigen::MatrixXd para, double posterior)
{
  if(posterior != numeric_limits<double>::infinity() 
     && posterior != -1*numeric_limits<double>::infinity())
    {
      //std::cout << "para = " << para <<"\n";

      if(std::isnan(posterior) == 0.0)
	{
	  for(unsigned int p=0; p<nPara; p++)
	    {
	      unsigned int id = hp.at(p).find_binID(para(0,p), gridsize, priorsMax.at(p),0);
	      double num = (double) hp.at(p).get_freq(id);
	      if(num == 0.0)
		hp.at(p).addY(id,posterior);
	      else
		{
		  hp.at(p).replaceY(id, hp.at(p).getY(id)*(num/(num+1))+posterior/(num+1));
		}
	      hp.at(p).add_freq(id);
	    }
	}
      else
	{
	  std::cout << "In marginal::add()\n";
	  std::cout << "posterior = NAN when para = " << para << "\n"; 
	}
    }
  return;
}

Eigen::MatrixXd marginal::compute_listJointPoints(std::vector<unsigned int> ids_para2integrate)
{

  // REMOVE
  /*
  for(unsigned int i=0; i<points2Integrate.size(); i++)
    {
      for(unsigned int j=0; j<points2Integrate.at(i).size(); j++)
	std::cout << points2Integrate.at(i).at(j) <<" ";
      std::cout << "\n";
    }
  */


	unsigned int nIDs_para = ids_para2integrate.size();
	unsigned int nRow = 1; //pow(nJointPoints,nIDs_para);
	for(unsigned int i=0; i< nIDs_para; i++)
	  {
	    unsigned int paraID = ids_para2integrate.at(i);
	    nRow *= points2Integrate.at(paraID).size();
	  }
	unsigned int nCol = nIDs_para;
	Eigen::MatrixXd list_otherPara(nRow,nCol);

	// REMOVE
	//std::cout << "nRow = " << nRow << " and nCol = " << nCol <<"\n";

	if(nCol == 1)
	{
	  for(unsigned int i=0; i<nJointPoints; i++)
	    {
	      // std::cout << "i = " << i <<"\n";
	      unsigned int id = ids_para2integrate.at(0);
	      list_otherPara(i,0) = points2Integrate.at(id).at(i);
	      if(priorsMax.at(id) == 0.0)
		i = nJointPoints;
	    }
	}
	else
	{
		std::vector<unsigned int> sub_ids = ids_para2integrate;
		sub_ids.erase(sub_ids.begin());
		Eigen::MatrixXd subList = compute_listJointPoints(sub_ids);
		
		// REMOVE
		//std::cout << "subList = " << subList <<"\n";
		
		unsigned int sub_nr  = subList.rows();
		for(unsigned int i=0; i<nJointPoints; i++)
		{
		  unsigned int id = ids_para2integrate.at(0);
		  list_otherPara.block(i*sub_nr,1,sub_nr,nCol-1) = subList;
		  for(unsigned int k=0; k<sub_nr; k++)
		    list_otherPara(i*sub_nr+k,0) = points2Integrate.at(id).at(i);
		  if(priorsMax.at(id) == 0.0)
		    i = nJointPoints;
		}
	}

	//REMOVE
	//std::cout << list_otherPara <<"\n\n";
	
	return list_otherPara;
}

void marginal::shareJointPoints()
{
  unsigned int nrow = pow(nJointPoints, nPara - nZeroPara);
  double point = 0.0;
  for(unsigned int i=0; i<nPara; i++){
    for(unsigned int j=0; j<nrow; j++){
      MPI::COMM_WORLD.Bcast(&jointDensity(j,i), 1,  MPI::DOUBLE, 0);
    }
  }
  return;
}

// YC 5/20/2014
// new version - keep the points that are commonly used for all the parameters.
void marginal::computeJointPoints()
{

  unsigned int nrow = pow(nJointPoints,nPara-nZeroPara);// + nPara*(gridsize-nJointPoints)*pow(nJointPoints,nPara-1);
  //	jointDensity.resize(nrow,nPara+1);
  //	jointDensity.setZero();
  
  std::vector<unsigned int> ids_allPara;
  for(unsigned int i=0; i<nPara; i++)
    ids_allPara.push_back(i);
  
  Eigen::MatrixXd subMat_points2Integrate = compute_listJointPoints(ids_allPara);
  
  jointDensity.block(0,0,nrow,nPara) = subMat_points2Integrate;
   
  return;
}



void marginal::computeJointDensity_MPI_overSample(IM im, popTree* poptree, Chain coldCh, unsigned int nProcs, unsigned int crr_procID)
{
  // REMOVE
  if(crr_procID == 0)
    std::cout << "\n In marginal::computeJointDensity_MPI_overSample()\n";
  // poptree->assign_populationSize(0.5);

  unsigned int nrow = jointDensity.rows();
  for(unsigned int i=0; i<nrow; i++)
    {
      Eigen::MatrixXd listPara= jointDensity.row(i);
      poptree->replacePara(listPara);
      // poptree->replacePara(im.get_ancPop(),listPara);
      

      // REMOVE
      // std::cout << "computing partial\n";
      double jointDensity_local = 0;
      jointDensity_local = coldCh.compute_partialJointPosteriorDensity_MPI_overSample(poptree, im, crr_procID, nProcs);
      // std::cout << "partialDensity = " << jointDensity_local << " on proc " << crr_procID <<"\n";
      double jointDensity_global = 0;
#ifdef MPI_ENABLED
      MPI::COMM_WORLD.Barrier();
      MPI::COMM_WORLD.Reduce(&jointDensity_local, &jointDensity_global, 1, MPI_DOUBLE, MPI_SUM, 0);
#endif
      if(crr_procID == 0)
	{
	  jointDensity(i,nPara) = jointDensity_global;
	  
	  poptree->print_allPopTree();
	  std::cout << "jointDensity_global = " << jointDensity_global <<"\n";
	}
    }
  return;
}

void marginal::computeJointDensity_MPI(IM im, popTree* poptree, Chain coldCh, unsigned int nProc, unsigned int crr_procID)
{

  unsigned int nrow = jointDensity.rows();
  for(unsigned int i=0; i<nrow; i++)
    {
      int procID = i-(unsigned int) floor((i)/nProc)*nProc;
      if(crr_procID == (unsigned) procID)
	{
	  
	  Eigen::MatrixXd listPara= jointDensity.row(i);
	  // Get a new poptree with the new parameters
	  poptree->replacePara(listPara);
	  // poptree->replacePara(im.get_ancPop(),listPara);
	  jointDensity(i,nPara) = coldCh.compute_jointPosteriorDensity_MPI(poptree, im,crr_procID);
	  if(procID != 0)
	    {
	      MPI::COMM_WORLD.Send(&jointDensity(i,nPara), 1, MPI::DOUBLE, 0, i);
	    }
	}
      else if(crr_procID == 0)
	{	  
	  MPI::Status status;
	  MPI::COMM_WORLD.Recv(&jointDensity(i,nPara), 1, MPI::DOUBLE, procID, i, status);
	}
    }  
  return;
}



void marginal::computeJointDensity(IM im, popTree* poptree, Chain coldCh)
{
	unsigned int nrow = jointDensity.rows();
	for(unsigned int i=0; i<nrow; i++)
	{
		Eigen::MatrixXd listPara= jointDensity.row(i);
		// Get a new poptree with the new parameters
		poptree->replacePara(listPara);
		// poptree->replacePara(im.get_ancPop(),listPara);
		jointDensity(i,nPara) = coldCh.compute_jointPosteriorDensity(poptree, im);
	}

	// REMOVE
	//std::cout << "jointDensity is\n" << jointDensity <<"\n\n";

	return;
}


/***
 * YC 11/18/2014
 * Parallelized over MCMC sample
 * New version - following how Ima2 computes the marginal distributions and means etc.
 */
void marginal::computeMarginalDensity_MPI_overSample(IM im, popTree* poptree, Chain coldCh,unsigned int nProcs, unsigned int crrProcID)
{

  std::vector<unsigned int> ids_para;
  for(unsigned int i=0; i<nPara; i++)
    ids_para.push_back(i);
  
  unsigned int nPoints = jointDensity.rows();
  
  if(crrProcID == 0)
    {
      for(unsigned int i=0; i<nPoints; i++)
	{
	  for(unsigned int p=0; p<nPara; p++)
	    {
	      if(i==0)
		{
		  // set the number of y values as 'gridsize' and initialize the values as zero - YC 7/17/2014
		  if(priorsMax.at(p) == 0.0)			  
		    hp.at(p).initiateY(1);
		  else
		    hp.at(p).initiateY(gridsize);
		}
	      unsigned int id = jointDensity(i,p)/binWidth_small.at(p)-0.5;
	      double val = jointDensity(i,nPara);
	      hp.at(p).addY((unsigned) id,val);
	    }
	}     
    }
  jointDensity.resize(0,0);
  
  Eigen::MatrixXi id(1,nPara);
  Eigen::MatrixXi id2int(1,nPara);
  
  // Compute the total number of grid points to compute
  unsigned int total_numJobs = 0;
  for(unsigned int p=0; p<nPara; p++)
    {
      unsigned int actual_gridsize = gridsize;
      if(priorsMax.at(p) ==0.0)
	{
	  actual_gridsize = 0;
	}
      total_numJobs += (actual_gridsize - nJointPoints);
    }
  total_numJobs *= pow(nJointPoints,nPara-1);

  unsigned int count_jobs =0;
  for(unsigned int p=0; p<nPara; p++)
    {
      unsigned int actual_gridsize = gridsize;
      if(priorsMax.at(p) ==0.0)
	{
	  actual_gridsize = 0;
	}
      for(unsigned int i=0; i< actual_gridsize; i++)
	{
	  
	  // YC 7/17/2017
	  // Modified the following statement. I guess it has a bug.
	  
	  double pseudoID = hp.at(p).getX(i)*nJointPoints/priorsMax.at(p)-0.5;
	  // double pseudoID = hp.at(p).getX(i)/binWidth_large.at(p)-0.5;	
	  // binWidth_large = (priorMax - priorMin)/nJointPoints
	  
	  if(abs(floor(pseudoID)-pseudoID)> 1/pow(10,6) && abs(ceil(pseudoID)-pseudoID)> 1/pow(10,6) ) // need to compute joint density
	    {
	      unsigned int sub_nr = 0;
	      Eigen::MatrixXd listPara; // containing the parameter values to compute
	      
	      std::vector<unsigned int> ids_subPara;
	      for(unsigned int pp=0; pp<nPara; pp++)
		{
		  if(pp != p)
		    {
		      ids_subPara.push_back(pp);
		    }
		}
	      Eigen::MatrixXd subMat = compute_listJointPoints(ids_subPara);
	      sub_nr = subMat.rows(); // "sub_nr" is equal to "(nPara-1)^nJointPoints".
	      listPara.resize(sub_nr,nPara);
	      for(unsigned int k=0; k<sub_nr; k++)
		{
		  listPara(k,p) = hp.at(p).getX(i);			    
		}
	      if(p!=0)
		{
		  listPara.block(0,0,sub_nr,p) =subMat.block(0,0,sub_nr,p);
		}
	      if(p!=nPara-1)
		{
		  listPara.block(0,p+1,sub_nr,nPara-p-1) = subMat.block(0,p,sub_nr,nPara-p-1);
		}
	      
	      for(unsigned int k=0; k<sub_nr; k++)
		{ 
		  Eigen::MatrixXd parameters = listPara.row(k);
		  // Get a new poptree with the new parameters
		  poptree->replacePara(parameters);
		  // poptree->replacePara(im.get_ancPop(),parameters);
		  double val_local = coldCh.compute_partialJointPosteriorDensity_MPI_overSample(poptree, im, crrProcID, nProcs);
		  double val_global = 0;
#ifdef MPI_ENABLED
		  MPI::COMM_WORLD.Barrier();
		  MPI::COMM_WORLD.Reduce(&val_local, &val_global, 1, MPI_DOUBLE, MPI_SUM, 0);
#endif
		  if(crrProcID == 0)
		    {
		      hp.at(p).addY(i,val_global);

		      std::cout << "\n In computeMarginalDensity_MPI_overSample()\n";
		      poptree->print_allPopTree();
		      std::cout << "val_global = " << val_global << "\n";
		    }
		  count_jobs++;
		}
	    }
	}
    }
  
  return;
}

/***
 * YC 5/20/2014
 * New version - following how Ima2 computes the marginal distributions and means etc.
 */
void marginal::computeMarginalDensity_MPI(IM im, popTree* poptree, Chain coldCh,unsigned int nProcs, unsigned int crrProcID)
{

  std::vector<unsigned int> ids_para;
  for(unsigned int i=0; i<nPara; i++)
    ids_para.push_back(i);
  
  unsigned int nPoints = jointDensity.rows();
  
  if(crrProcID == 0)
    {
      for(unsigned int i=0; i<nPoints; i++)
	{
	  for(unsigned int p=0; p<nPara; p++)
	    {
	      if(i==0)
		{
		  // set the number of y values as 'gridsize' and initialize the values as zero - YC 7/17/2014
		  if(priorsMax.at(p) == 0.0)			  
		    hp.at(p).initiateY(1);
		  else
		    hp.at(p).initiateY(gridsize);
		}
	      unsigned int id = jointDensity(i,p)/binWidth_small.at(p)-0.5;
	      double val = jointDensity(i,nPara);
	      hp.at(p).addY((unsigned) id,val);
	    }
	}      
     
    }
  jointDensity.resize(0,0);
  
  Eigen::MatrixXi id(1,nPara);
  Eigen::MatrixXi id2int(1,nPara);
  //double prod_binWidthLarge = 1.0;
  
  // Compute the total number of grid points to compute
  unsigned int total_numJobs = 0;
  for(unsigned int p=0; p<nPara; p++)
    {
      unsigned int actual_gridsize = gridsize;
      if(priorsMax.at(p) ==0.0)
	{
	  actual_gridsize = 0;
	}
      total_numJobs += (actual_gridsize - nJointPoints);
    }
  total_numJobs *= pow(nJointPoints,nPara-1);

  unsigned int count_jobs =0;
  for(unsigned int p=0; p<nPara; p++)
    {
      unsigned int actual_gridsize = gridsize;
      if(priorsMax.at(p) ==0.0)
	{
	  actual_gridsize = 0;
	}
      for(unsigned int i=0; i< actual_gridsize; i++)
	{
	  
	  // YC 7/17/2017
	  // Modified the following statement. I guess it has a bug.
	  
	  double pseudoID = hp.at(p).getX(i)*nJointPoints/priorsMax.at(p)-0.5;
	  // double pseudoID = hp.at(p).getX(i)/binWidth_large.at(p)-0.5;	
	  // binWidth_large = (priorMax - priorMin)/nJointPoints
	  
	  // REMOVE
	  //std:: cout << "pseudoID = " << pseudoID << " floor = "<< floor(pseudoID) << " ceiling = " << ceil(pseudoID) <<"\n";
	  
	  if(abs(floor(pseudoID)-pseudoID)> 1/pow(10,6) && abs(ceil(pseudoID)-pseudoID)> 1/pow(10,6) ) // need to compute joint density
	    {
	      unsigned int sub_nr = 0;
	      Eigen::MatrixXd listPara; // containing the parameter values to compute
	      
	      std::vector<unsigned int> ids_subPara;
	      for(unsigned int pp=0; pp<nPara; pp++)
		{
		  if(pp != p)
		    {
		      ids_subPara.push_back(pp);
		    }
		}
	      Eigen::MatrixXd subMat = compute_listJointPoints(ids_subPara);
	      sub_nr = subMat.rows(); // "sub_nr" is equal to "(nPara-1)^nJointPoints".
	      listPara.resize(sub_nr,nPara);
	      for(unsigned int k=0; k<sub_nr; k++)
		{
		  listPara(k,p) = hp.at(p).getX(i);			    
		}
	      if(p!=0)
		{
		  listPara.block(0,0,sub_nr,p) =subMat.block(0,0,sub_nr,p);
		}
	      if(p!=nPara-1)
		      {
			listPara.block(0,p+1,sub_nr,nPara-p-1) = subMat.block(0,p,sub_nr,nPara-p-1);
		      }
		  	    
		    for(unsigned int k=0; k<sub_nr; k++)
		      { 
			
			int procID = count_jobs - nProcs*floor(count_jobs/nProcs);
			// int procID = i*nPara+p - floor((i*nPara+p)/nProcs)*nProcs;

		
			double val = 0.0;
			if(procID == crrProcID)
			  {
			    Eigen::MatrixXd parameters = listPara.row(k);
			    // Get a new poptree with the new parameters
			    poptree->replacePara(parameters);
			    // poptree->replacePara(im.get_ancPop(),parameters);
			    val = coldCh.compute_jointPosteriorDensity_MPI(poptree, im,  crrProcID);
			 
			    if(crrProcID != 0)
			      {
				//std::cout << "send val = "<< val <<" with tag= "<< k <<" from proc "<< procID <<" p= "<< p <<", i="<< i <<"\n";
				MPI::COMM_WORLD.Send(&val, 1, MPI::DOUBLE, 0, k);				  
			      }
			    else
			      {
				hp.at(p).addY(i,val);
			      }
			  }
			if(crrProcID == 0 && crrProcID != procID)
			  {
			    //std::cout << "Will receive val with tag= "<< k <<" from proc "<< procID <<" p= "<< p <<", i="<< i <<"\n";
			    MPI::Status status;
			    MPI::COMM_WORLD.Recv(&val, 1, MPI::DOUBLE, procID, k, status);
			    //std::cout << "Received val = "<< val <<" with tag= "<< k <<" from proc "<< procID <<" p= "<< p <<", i="<< i <<"\n";
			    hp.at(p).addY(i,val);
			  }	
			count_jobs++;
		      }
		  }
	      }
	  }
	
	return;
}


/***
 * YC 5/20/2014
 * New version - following how Ima2 computes the marginal distributions and means etc.
 */
void marginal::computeMarginalDensity(IM im, popTree* poptree, Chain coldCh)
{
	std::vector<unsigned int> ids_para;
	for(unsigned int i=0; i<nPara; i++)
		ids_para.push_back(i);

	unsigned int nPoints = jointDensity.rows();


	for(unsigned int i=0; i<nPoints; i++)
	{
		for(unsigned int p=0; p<nPara; p++)
		{
			if(i==0)
				hp.at(p).initiateY(gridsize);
			unsigned int id = jointDensity(i,p)/binWidth_small.at(p)-0.5;
			double val = jointDensity(i,nPara);
			hp.at(p).addY((unsigned) id,val);
		}
	}
	jointDensity.resize(0,0);

	Eigen::MatrixXi id(1,nPara);
	Eigen::MatrixXi id2int(1,nPara);
	double pseudoID = 0.0;
	//double prod_binWidthLarge = 1.0;
	for(unsigned int i=0; i<gridsize; i++)
	{
		for(unsigned int p=0; p<nPara; p++)
		{
			pseudoID = hp.at(p).getX(i)/binWidth_large.at(p)-0.5;
			if(floor(pseudoID) != ceil(pseudoID)) // need to compute joint density
			{
				std::vector<unsigned int> ids_subPara;
				for(unsigned int pp=0; pp<nPara; pp++)
					if(pp != p)
						ids_subPara.push_back(pp);
				Eigen::MatrixXd subMat = compute_listJointPoints(ids_subPara);
				unsigned int sub_nr = subMat.rows();
				Eigen::MatrixXd listPara(sub_nr,nPara);
				for(unsigned int k=0; k<sub_nr; k++)
					listPara(k,p) = hp.at(p).getX(i);
				if(p!=0)
					listPara.block(0,0,sub_nr,p) =subMat.block(0,0,sub_nr,p);
				if(p!=nPara-1)
					listPara.block(0,p+1,sub_nr,nPara-p-1) = subMat.block(0,p,sub_nr,nPara-p-1);


				for(unsigned int k=0; k<sub_nr; k++)
				{
					Eigen::MatrixXd parameters = listPara.row(k);
					// Get a new poptree with the new parameters
					// poptree->replacePara(im.get_ancPop(),parameters);
					poptree->replacePara(parameters);
					double val = coldCh.compute_jointPosteriorDensity(poptree, im);
					hp.at(p).addY(i,val);
				}
			}
		}
	}

	return;
}

void hist::compute_bin_MeanPosterior_forDE()
{
  // std::cout << "In hist::compute_bin_MeanPosterior_forDE()\n";
  //std::cout << "y=" << y <<"\n";
  //std::cout << "freq=" << freq <<"\n\n";

  for(unsigned int i=0; i<y.cols(); i++)
    {
      if(freq(0,i) !=0)
	y(0,i) = y(0,i)/freq(0,i);
    }

  return;
}

void hist::computeMeanIntervals()
{
  /*
  std::cout << "in hist::computeMeanIntervals()\n";
  std::cout << "y = " << y <<"\n";
  std::cout << "x = " << x <<"\n";
  std::cout << "freq = " << freq <<"\n";
  */

  double sumY = y.sum();
  double unif = 1/(double) x.cols();
	//p = y/y.sum();
  Eigen::MatrixXd xy = x.cwiseProduct(y);
	//Eigen::MatrixXd xp = x.cwiseProduct(p); // xy.sum()/sumY;
  if(sumY == 0)
    {	    
      mean = x.sum()* unif;
    }
  else
    mean = xy.sum()/sumY;


  unsigned int count =0;
  unsigned int found =0;
  double cumProb =0.0;
  while(count < y.cols() && found == 0)
    {
      if(sumY == 0)
	cumProb += unif;
      else
	cumProb += y(0,count)/sumY;
      
      if(cumProb >= 0.025)
	{
	  found = 1;
	  lower95 = x(0,count);
	}
      count++;
    }
  
  count = y.cols()-1;
  found =0;
  cumProb = 0;
  while(count>=0 && found ==0)
    {
      if(sumY == 0)
	cumProb += unif;
      else
	cumProb += y(0,count)/sumY;
      if(cumProb >=0.025)
	{
	  found = 1;
	  upper95 = x(0,count);
	}
      count--;
    }
  return;
}

void marginal::computeMeanIntervals_forDE()
{
  //std::cout << "In marginal::computeMeanIntervals_forDE()\n";
  for(unsigned int p=0; p<nPara; p++)
    {
      //std::cout << "for p = " << p <<" computing mean posterior of each bin\n";
      //hp.at(p).compute_bin_MeanPosterior_forDE();
      //std::cout << "computing the posterior means\n";
      hp.at(p).computeMeanIntervals();
    }
  return;
}

void marginal::computeMeanIntervals()
{
  for(unsigned int p=0; p<nPara; p++)
    {
      hp.at(p).computeMeanIntervals();
    }
  return;
}

void marginal::printMeanIntervals()
{
	for(unsigned int i=0; i<nPara; i++)
	{
	  /*
		if(i<3)
			std::cout << "Population size N" << i+1 <<":\n";
		else if(i<5)
			std::cout << "Migration rate m" << i-2 <<":\n";
		else
			std::cout <<"Splitting time:\n";
	  */
	  std::cout << "\tMean: " << hp.at(i).getMean()<<"\n";
	  std::cout << "\t95% interval (" << hp.at(i).getLower95() <<"," << hp.at(i).getUpper95() <<")\n\n";
	}
}


void marginal::saveMeanIntervals()
{
  ofstream file;
  file.open ("mean_intervals.txt");
  for(unsigned int i=0; i<nPara; i++)
    {
      if(i<3)
	file << "Population size N" << i+1 <<":\n";
      else if(i<5)
	file << "Migration rate m" << i-2 <<":\n";
      else
	file <<"Splitting time:\n";
      file << "\tMean: " << hp.at(i).getMean()<<"\n";
      file << "\t95% interval (" << hp.at(i).getLower95() <<"," << hp.at(i).getUpper95() <<")\n";
    }
  file.close();
  return;
}

void marginal::printHistogram()
{ 
  
  for(unsigned int p =0;p<6; p++)
    {
      std::cout << "Histogram for p = "<< p <<"\n";
      std::cout << "x = "<< hp.at(p).getX() <<"\n" <<"y=" << hp.at(p).getY() <<"\n";
    }
  std::cout << "\n";
  for(unsigned int p =0;p<6; p++)
    {
      std::cout << "Histogram (prob) for p = "<< p <<"\n";
      Eigen::MatrixXd y_p = hp.at(p).getY();
      if(y_p.sum() == 0)
	{
	  double unif = 1/(double) gridsize;
	  std::cout  << "x = "<< hp.at(p).getX() <<"\n" <<"p=";
	  for(unsigned int i=0; i<gridsize; i++)
	    std::cout << unif <<" ";
	  std::cout <<"\n";
	}
      else
	std::cout << "x = "<< hp.at(p).getX() <<"\n" <<"p=" << y_p/y_p.sum() <<"\n";
    }
 return;
}

void marginal::saveHistogram_forDE()
{ 
  ofstream file;
  file.open ("poptree_Histogram.txt");

  for(unsigned int p=0; p<nPara; p++)
    {
      file << "p" << p+1 << " prob freq "; 
    }
  file <<"\n";

  for(unsigned int i=0; i<gridsize; i++)
    {
      for(unsigned int p=0; p<nPara; p++)
	{
	  file << hp.at(p).getX(i) << " " 
	       << hp.at(p).getY(i) << " "
	       << hp.at(p).get_freq(i) << " ";
	}
      file <<"\n";
    }
  /*
  for(unsigned int p =0;p<6; p++)
    {
      file << "Histogram for p = "<< p <<"\n";
      // hp.at(p).printHistorgram(); 
      file << "x = "<< hp.at(p).getX() <<"\n" <<"y=" << hp.at(p).getY() <<"\n";
    }
  file << "\n";
  for(unsigned int p =0;p<6; p++)
    {
      file << "Histogram (prob) for p = "<< p <<"\n";
      Eigen::MatrixXd y_p = hp.at(p).getY();
      if(y_p.sum() == 0)
	{
	  double unif = 1/(double) gridsize;
	  file << "x = "<< hp.at(p).getX() <<"\n" <<"p=";
	  for(unsigned int i=0; i<gridsize; i++)
	    file << unif <<" ";
	  file <<"\n";
	}
      else
	file << "x = "<< hp.at(p).getX() <<"\n" <<"p=" << y_p/y_p.sum() <<"\n";
    }
  */
  file.close();
 return;
}

void marginal::saveHistogram()
{ 
  ofstream file;
  file.open ("poptree_Histogram.txt");
  
  for(unsigned int p =0;p<6; p++)
    {
      file << "Histogram for p = "<< p <<"\n";
      // hp.at(p).printHistorgram(); 
      file << "x = "<< hp.at(p).getX() <<"\n" <<"y=" << hp.at(p).getY() <<"\n";
    }
  file << "\n";
  for(unsigned int p =0;p<6; p++)
    {
      file << "Histogram (prob) for p = "<< p <<"\n";
      Eigen::MatrixXd y_p = hp.at(p).getY();
      if(y_p.sum() == 0)
	{
	  double unif = 1/(double) gridsize;
	  file << "x = "<< hp.at(p).getX() <<"\n" <<"p=";
	  for(unsigned int i=0; i<gridsize; i++)
	    file << unif <<" ";
	  file <<"\n";
	}
      else
	file << "x = "<< hp.at(p).getX() <<"\n" <<"p=" << y_p/y_p.sum() <<"\n";
    }
  file.close();
 return;
}

void marginal::prepareComputingHistrograms_allProc(IM im)
{
  initiate(im);
  initiate_histXvalues_binWidth_small();
  initiate_Points2Integrate_largeBinWidth();
}

void marginal::prepareComputingHistrograms(IM im)
{
  computeJointPoints();
  //std::cout << "Done2 - remove it later\n";
}
