#include "UGRNN.h"

UGRNN::UGRNN(Molecule** moleculeList,myOptions & _Opt ,int moleculeNumber,string modelDir,int predict){
	setModelDir(modelDir);
	setPropertiesNum(_Opt.getProperties());
	if(predict==0){
		setLogTP(modelDir+"/logTP");
		setLogTPVal(modelDir+"/logTPVal");
		totalAtoms = getAllAtoms(moleculeList,moleculeNumber);
		averageAtomsNumber = computeAverageAtomsNumber(moleculeList,moleculeNumber);
		writeTotalAtoms(modelDir+"/totalAtoms");
		Opt = _Opt;
		cout << "Seed=" << Opt.getSeed() << endl;
		//@debug
	    for(int i=0;i<totalAtoms.size();i++){
			cout << totalAtoms.at(i) << " ";
		}
		cout << endl;
		///////
		ATOM_NUMBER = totalAtoms.size();
		MAX_CHILDREN = getMaxNumberOfChildren(moleculeList,moleculeNumber);
		writeMAX_CHILDREN(modelDir+"/MAX_CHILDREN");
		//debug
		cout << "MAX_CHILDREN=" << MAX_CHILDREN << endl;	
		///////////////
	}else{
		vector<string> totalAtomsPredict;
		int MAX_CHILDRENPredict;
		totalAtomsPredict = getAllAtoms(moleculeList,moleculeNumber);
		readTotalAtoms(modelDir+"/totalAtoms");
		compareTotalAtoms(totalAtoms,totalAtomsPredict);
		ATOM_NUMBER = totalAtoms.size();
		MAX_CHILDRENPredict = getMaxNumberOfChildren(moleculeList,moleculeNumber);
		readMAX_CHILDREN(modelDir+"/MAX_CHILDREN");
		if(MAX_CHILDREN < MAX_CHILDRENPredict){
			cout << "MAX_CHILDREN=" << MAX_CHILDREN << endl;	
			cout << "MAX_CHILDRENPredict=" << MAX_CHILDRENPredict << endl;
			cout << "totalAtoms=" << totalAtoms.size() << endl;
			cout << "totalAtomsPredict=" << totalAtomsPredict.size() << endl;	
			exit(1);		
		}
		Opt = _Opt;
		//debug
		cout << "MAX_CHILDREN=" << MAX_CHILDREN << endl;	
		//cout << "MAX_CHILDRENPredict:" << MAX_CHILDRENPredict << endl;	
		///////////////
	}
	//1 atom label + 1 bond type + G_dimension * output_number 
	TOTAL_INPUTS = 1+MAX_CHILDREN*Opt.getMNN_OUTPUTS_NUMBER()+MAX_CHILDREN;
	cout << "TOTAL_INPUTS=" << TOTAL_INPUTS << endl;
	cout << "ATOM_NUMBER=" << ATOM_NUMBER << endl; 
	//initialize MNN
	//total input:1 categorical + MAX_CHILDREN*getMNN_OUTPUTS_NUMBER() real valued
	//MNN isn't an output network
	int* MNNk = new int[TOTAL_INPUTS];
	MNNk[0] = ATOM_NUMBER;
	for(int i=1;i<TOTAL_INPUTS;i++){
		MNNk[i]=1;
		
	}
	MNN = new NNt(1,(MAX_CHILDREN*Opt.getMNN_OUTPUTS_NUMBER())+MAX_CHILDREN,Opt.getMNN_HIDDEN_UNITS(),Opt.getMNN_OUTPUTS_NUMBER(),MNNk,0,1);			
	MNN->initWeights(Opt.getSeed());
	//initialize ONN 
	int* ONNk;
	if(propertiesNum<=0){
		ONNk = new int[Opt.getMNN_OUTPUTS_NUMBER()];	
		for(int i=0;i<(Opt.getMNN_OUTPUTS_NUMBER());i++){
			ONNk[i] = 1;
		}
		ONN = new NNt(0,Opt.getMNN_OUTPUTS_NUMBER(),Opt.getONN_HIDDEN_UNITS(),Opt.getONN_OUTPUTS_NUMBER(),ONNk,1,1); 
		ONN->initWeights(Opt.getSeed());
	}else if(propertiesNum == 1){
		cout << "popertiesNum=" << propertiesNum << endl;
		ONNk = new int[Opt.getMNN_OUTPUTS_NUMBER()+1];	
		ONNk[0] = propertiesNum;
		//ONNk[0]=0;
		for(int i=1;i<(Opt.getMNN_OUTPUTS_NUMBER()+1);i++){
			ONNk[i] = 1;
		}
		ONN = new NNt(1,Opt.getMNN_OUTPUTS_NUMBER(),Opt.getONN_HIDDEN_UNITS(),Opt.getONN_OUTPUTS_NUMBER(),ONNk,1,1); 
		ONN->initWeights(Opt.getSeed());
		
	}else{
		exit(1);
	}
	//ONN = new NNt(1,Opt.getMNN_OUTPUTS_NUMBER(),Opt.getONN_HIDDEN_UNITS(),Opt.getONN_OUTPUTS_NUMBER(),ONNk,1,1); 
	delete[] MNNk;
	delete[] ONNk;
	
}

UGRNN::~UGRNN(){
	cout << "UGRNN DESTRUCTOR" << endl;
	delete ONN;
	delete MNN;
}

void UGRNN::write(ofstream & model,ofstream & modelM){
	ONN->write(model); 	
	model.close();
	MNN->write(modelM);
	modelM.close();
}

void UGRNN::read(string ONNModelPath,string MNNModelPath){
	ifstream ONNModel,MNNModel;
	ONNModel.open(ONNModelPath.c_str());
	MNNModel.open(MNNModelPath.c_str());
	if(!ONNModel.is_open()){
		cerr << "CANNOT READ ONN MODEL" << endl;
		exit(1);
	}
	if(!MNNModel.is_open()){
		cerr << "CANNOT OPEN MNN MODEL" << endl;
		exit(1);
	}
	ONN->read(ONNModel);
	MNN->read(MNNModel);
	MNNModel.close();	
	ONNModel.close();
}

Float UGRNN::forward(Molecule* molecule,Float* properties,Float* descriptors){
	//order graph and load neighbours
	molecule->g->orderGraph();
	molecule->g->loadNeighbours();
	/////////////////////////////////	
	//cout << "Molecule " << i << endl;
	Float* ONNOut = new Float[1];
	memset(ONNOut,0,sizeof(Float));
	Float* MNNOut = new Float[Opt.getMNN_OUTPUTS_NUMBER()];
	memset(MNNOut,0,sizeof(Float)*Opt.getMNN_OUTPUTS_NUMBER());
	//forward MNN
	NodeIO*** nodeIO = new NodeIO**[molecule->g->total_vert];
	initNodeIO(molecule,nodeIO);
	forwardMNN(molecule,nodeIO,MNNOut);		
	/////////////
	//forward ONN
	//debug
	//cout << "properties[0]=" << properties[0] << endl; 
	///////	
	ONN->forward(properties,MNNOut);
	//print MNN output
	//if(e==Opt.getNEpochs()-1){
	//	for(int d=0;d<(MNN->get_NO());d++){
	//		logTestVector << MNN->out()[d] << " ";
	//	}
	//	logTestVector << t[i];
	//	if(i<(validationVector.size())-1){
	//		logTestVector << endl;
	//	}
	//}
	//////////////////
	/////////////
	//assing MNNOut to descriptors
	if(descriptors!=NULL){	
		for(int d=0;d<MNN->get_NO();d++){
			descriptors[d]=MNNOut[d];
		}
	}
	//////////////////////////////
	//free memory
	delete [] MNNOut;
	delete [] ONNOut;
	for(int a=0;a<molecule->g->total_vert;a++){
		for(int b=0;b<molecule->g->total_vert;b++){
			delete nodeIO[a][b];
		}
		delete [] nodeIO[a];
	}
	delete [] nodeIO; 
	/////////////
	//unload graph
	molecule->g->unloadGraph();	
	//////////////

	return ONN->out()[0];
}

void UGRNN::forwardMNN(Molecule* molecule,NodeIO*** nodeIO,Float* out){
//iterate for each subgraph j in the molecule
for(int j=0;j<molecule->g->total_vert;j++){
//for(int j=0;j<1;j++){
	//@debug
	//cout << "forwarding  graph: " << j << endl;
	for(int k=0;k<molecule->g->total_vert;k++){
		//set I2 elements to 0	
		nodeIO[j][k]->setI2();
	}
	////////
	//@debug
	//for(int k=0;k<moleculeList[molPosition[i]]->g->total_vert;k++){
	//	nodeIO[j][k]->printI1();
	//}
	////////
	sum(out,compute(j,nodeIO[j],molecule->g->children[j],molecule->g->num_children[j],molecule->g->num_children[j][j],molecule->getBondType()),MNN->get_NO());
	//normalizeOut(out,MNN->get_NO(),10);
}
/*	
for(int j=0;j<MNN->get_NO();j++){
	out[j]=tanh(out[j]);
}
*/
	
}

void UGRNN::initNodeIO(Molecule* molecule,NodeIO*** nodeIO){
	for(int a=0;a<molecule->g->total_vert;a++){
		nodeIO[a]=new NodeIO*[molecule->g->total_vert];
		for(int b=0;b<molecule->g->total_vert;b++){
			nodeIO[a][b] = new NodeIO(ATOM_NUMBER,MAX_CHILDREN*MNN->get_NO()+MAX_CHILDREN,MNN->get_NO());
			//@debug	
			//cout << moleculeList[molPosition[i]]->g->vert[a]->symbol << " mapped as: " <<  unTranslateU(moleculeList[molPosition[i]]->g->vert[a]->symbol) << endl;
			////////
			nodeIO[a][b]->categoricalInit(unTranslateU(molecule->g->vert[b]->symbol),totalAtoms);
			//@debug
			//nodeIO[a][b]->printI1();
			////////	
		}
	}
}

void UGRNN::backwardMNN(Molecule* molecule,NodeIO*** nodeIO){
	//backprop for each subgraph
	Float* delta = new Float[Opt.getMNN_OUTPUTS_NUMBER()];
	//cout << "delta[0] = " << ONN->back_out()[0] << endl;
	for(int k=propertiesNum;k<(Opt.getMNN_OUTPUTS_NUMBER()+propertiesNum);k++){
		delta[k-propertiesNum] = ONN->back_out()[k];
		//debug
		//cout << delta[k] << " ";
		///////	
	}	
	//debug
	//cout << endl;
	///////
	for(int j=0;j<molecule->g->total_vert;j++){
	//for(int j=0;j<1;j++){
		computeBack(j,nodeIO[j],delta,molecule->g->children[j],molecule->g->num_children[j],molecule->g->num_children[j][j]);
		//MNN->scaleGradient((Float)molecule->g->total_vert);
	}
	delete[] delta;
}


void UGRNN::trainAndValidate(Molecule** train,Molecule** validation,const vector<Float> & bound,int numberOfModels){
	//global variables
	int shuffleSeed = 5;
	vector<int> molPosition;
	Float* dummy = new Float[1];
	memset(dummy,0,sizeof(Float));
	Float epsilon = Opt.getEpsilon();
	Float epsilonMin = Opt.getEpsilonMin();
	Float epStep = (epsilon-epsilonMin)/(Float)Opt.getEpsEpoch();
	Float* lowestAAE = new Float[numberOfModels];
	memset(lowestAAE,0,sizeof(Float)*numberOfModels);
	Float** AAE = new Float*[Opt.getNEpochs()];	
	for(int i=0;i<Opt.getNEpochs();i++){
		AAE[i] = new Float[2];
	}
	Float CPU_time_usage;
	//initialize predicted vector and target vector
	vector<Float> bestP(validationVector.size(),0);
	vector<Float> bestT(validationVector.size(),0);
	///////////////////////////////////////////////
	//////////////////
	//initialize epsilon
    for(int e=0;e<firstEpoch;e++){
        epsilon-=epStep;
    }
    ////////////////////
    //////////////////
    for(int e=firstEpoch;e<Opt.getNEpochs();e++){
		cout << "Epoch:" << e << endl;
		//cout << "lowestAAE: " << lowestAAE << endl;
		molPosition = shuffle(trainVector.size(),shuffleSeed);	
		//debug
		//vector<int>::iterator it = molPosition.begin();
		//while(it!=molPosition.end()){
		//	cout << *it++ << " ";
		//}
		//cout << endl;
		///////
		//initialize predicted vector and target vector
		vector<Float> p(trainVector.size(),0);
		vector<Float> t(trainVector.size(),0); 
		///////////////////////////////////////////////
		//init error vector and pearson variable
		vector<Float> error(2,0);
		Float PEARSON = 0;
		////////////////////////////////////////
		clock_t epoch_begin = clock();
		for(int i=0;i<trainVector.size();i++){			
			//cout << "Molecule " << i << endl; 			
			//order graph and load neighbours			
			clock_t begin = clock();			
			train[molPosition[i]]->g->orderGraph();	
			train[molPosition[i]]->g->loadNeighbours();
			clock_t end = clock();			
			//CPU_time_usage = get_CPU_time_usage(end,begin);
			//cout << "CPU_time_usage for loading graph = " << CPU_time_usage << endl;
			//cout << "Graph loaded" << endl;
			/////////////////////////////////
			begin = clock();			
			Float* ONNOut = new Float[1];
			memset(ONNOut,0,sizeof(Float));
			Float* target = new Float[1];
			memset(target,0,sizeof(Float));
			Float* MNNOut = new Float[Opt.getMNN_OUTPUTS_NUMBER()];
			memset(MNNOut,0,sizeof(Float)*Opt.getMNN_OUTPUTS_NUMBER());
			//forward through MNN
			NodeIO*** nodeIO = new NodeIO**[train[molPosition[i]]->g->total_vert];
			initNodeIO(train[molPosition[i]],nodeIO); 
			forwardMNN(train[molPosition[i]],nodeIO,MNNOut);
			/////////////////////	
			//forward through ONN
			dummy[0] = (Float)train[molPosition[i]]->logP;
			ONN->forward(dummy,MNNOut);
			ONNOut[0] = ONN->out()[0];
			/////////////////////
			//init target
			target[0] = train[molPosition[i]]->target;
			//target[0] = train[molPosition[i]]->trueTarget;
			/////////////
			//backward through ONN
			ONN->backward(target);
			ONN->gradient(dummy,MNNOut,target);
			//////////////////////
			//backward through MNN
			backwardMNN(train[molPosition[i]],nodeIO);
			//////////////////////
			//save p and t values
			p[i] = denormalize(ONN->out()[0],bound);
			//p[i] = ONN->out()[0];
			//t[i] = denormalize(target[0],bound);
			t[i] = train[molPosition[i]]->trueTarget;
			/////////////////////
			//update weights
			if(i%Opt.getBatchBlocks()==0 && i!=0){
				ONN->updateWeightsClipped(epsilon);
				MNN->updateWeightsClipped(epsilon);
				ONN->resetGradient();
				MNN->resetGradient();
			}//else if(i%Opt.getBatchBlocks()==0 && i!=0 && e>1000){
			//	if (e==1001){
			//		ONN->initWeights(Opt.getSeed());	
			//	} 
			//	ONN->updateWeightsClipped(epsilon);
			//	ONN->resetGradient();
			//}	
			////////////////
			//free memory
			delete[] ONNOut;
			delete[] target;
			delete[] MNNOut;
			for(int a=0;a<train[molPosition[i]]->g->total_vert;a++){
				for(int b=0;b<train[molPosition[i]]->g->total_vert;b++){
					delete nodeIO[a][b];
				}
				delete[] nodeIO[a];
			}
			delete[] nodeIO;	
			/////////////
			//unload graph
			train[molPosition[i]]->g->unloadGraph();
			//////////////
			end = clock();
			//CPU_time_usage = get_CPU_time_usage(end,begin);
			//cout << "molecule CPU_time_usage = " << CPU_time_usage << endl;
		}
		//update weights for last molecules
		//if(e<=1000){
			ONN->updateWeightsClipped(epsilon);
			MNN->updateWeightsClipped(epsilon);
			ONN->resetGradient();
			MNN->resetGradient();
		//}else if(e>1000){
		//	ONN->updateWeightsClipped(epsilon);
		//	ONN->resetGradient();
		//}
		///////////////////////////////////
		//display error
		//setLogTP(modelDir+"/logTP");
		setLogError(modelDir+"/logError");
		error = computeError(t,p,trainVector.size());
		PEARSON = computeCorrelation(t,p,trainVector.size());
		cout << "RMSE and AAE at epoch " << e << " " << error[0] << " " << error[1] << endl;
		cout << "PEARSON " << PEARSON << endl;
		cout << "epsilon:" << epsilon << endl;
		logError << error[0] << "\t" << error[1] << "\t" << PEARSON <<  endl; 
		///////////////
		//update epsilon and shuffleSeed
		epsilon -= epStep;
		shuffleSeed++;
		////////////////////////////////
		//validate
		this->validate(validation,bound,lowestAAE,bestP,bestT,e,numberOfModels,AAE);	
		//////
		//close logFiles
		logError.close();
		////////////////
		clock_t epoch_end = clock();
		CPU_time_usage = get_CPU_time_usage(epoch_end,epoch_begin);
		cout << "epoch CPU_time_usage = " << CPU_time_usage << endl;	
	}
	//free memory
	delete[] dummy;
	delete[] lowestAAE;
	for(int i=0;i<Opt.getNEpochs();i++){
		delete[] AAE[i];
	}
	delete[] AAE;
	/////////////	
}

void UGRNN::validate(Molecule** validation,const vector<Float> & bound,Float* lowestAAE,vector<Float> & bestP,vector<Float> & bestT,int e,int numberOfModels,Float** AAE){
	//global variables
	Float* dummy = new Float[4];
	memset(dummy,0.0,sizeof(Float)*4);
	vector<Float> error(2,0);
	Float PEARSON = 0;
        //////////////////
	setLogValidationError(modelDir+"/logValidationError");
	//read model
	//debug
	//cout << modelDir+"/ONN" << endl;
	//cout << modelDir+"/MNN" << endl;
	///////
	//MNN->initWeights(Opt.getSeed());
	//ONN->initWeights(Opt.getSeed());
	//read(modelDir+"/ONN",modelDir+"/MNN");	
	////////////
	//init target and predicted vectors
	vector<Float> t(validationVector.size(),0);
	vector<Float> p(validationVector.size(),0);
	///////////////////////////////////
	for(int i=0;i<validationVector.size();i++){
		//order graph and load neighbours			
		validation[i]->g->orderGraph();	
		validation[i]->g->loadNeighbours();
		/////////////////////////////////
		//cout << "Molecule " << i << endl;
		Float* ONNOut = new Float[1];
		memset(ONNOut,0,sizeof(Float));
		Float* target = new Float[1];
		memset(target,0,sizeof(Float));
		Float* MNNOut = new Float[Opt.getMNN_OUTPUTS_NUMBER()];
		memset(MNNOut,0,sizeof(Float)*Opt.getMNN_OUTPUTS_NUMBER());
		//forward MNN
		NodeIO*** nodeIO = new NodeIO**[validation[i]->g->total_vert];
		initNodeIO(validation[i],nodeIO);
		forwardMNN(validation[i],nodeIO,MNNOut);		
		/////////////
		//forward ONN
		dummy[0] = (Float)validation[i]->logP;
		ONN->forward(dummy,MNNOut);
		ONNOut[0] = ONN->out()[0];
		/////////////
		//init target
		target[0] = validation[i]->target;
		//target[0] = validation[i]->trueTarget;
		/////////////
		//save target and predicted values
		p[i] = denormalize(ONN->out()[0],bound);
		t[i] = validation[i]->trueTarget;
		//////////////////////////////////
		//print MNN output
		if(e==Opt.getNEpochs()-1){
			for(int d=0;d<(MNN->get_NO());d++){
				logTestVector << MNN->out()[d] << " ";
			}
			logTestVector << t[i];
			if(i<(validationVector.size())-1){
				logTestVector << endl;
			}
		}
		//////////////////
		//free memory
		delete [] MNNOut;
		delete [] ONNOut;
		delete [] target;
		for(int a=0;a<validation[i]->g->total_vert;a++){
			for(int b=0;b<validation[i]->g->total_vert;b++){
				delete nodeIO[a][b];
			}
			delete [] nodeIO[a];
		}
		delete [] nodeIO; 
		/////////////
		//unload graph
		validation[i]->g->unloadGraph();
		//////////////	
	}
	//display error
	error = computeError(t,p,validationVector.size());
	if(validationVector.size()>1){
		PEARSON = computeCorrelation(t,p,validationVector.size());
		cout << "Validation RMSE and AAE at epoch " << e << " " << error[0] << " " << error[1] << endl;
		cout << "Validation PEARSON " << PEARSON << endl;
		//update AAE
		AAE[e][0] = e;
		AAE[e][1] = error[0];
		////////////
		//save model
		ofstream MNNModel, ONNModel;
		//es stands for epoch stream 
		stringstream es;
		es << e;
		MNNModel.open((modelDir+"/MNN"+es.str()).c_str());
		ONNModel.open((modelDir+"/ONN"+es.str()).c_str());
		if(!MNNModel.is_open()){
			cerr << "CANNOT WRITE MNN MODEL" << endl;	
		}
		if(!ONNModel.is_open()){
			cerr << "CANNOT WRITE ONN MODEL" << endl;
		}
		write(ONNModel,MNNModel);	
		////////////
		logValidationError << error[0] << "\t" << error[1] << "\t" << PEARSON <<  endl; 
	}else{
		//this part is unusefull because validationVector.size() cannot be = 1
		cout << "Validation RMSE and AAE at epoch " << e << " " << error[0] << " " << error[1] << endl;
		logValidationError << error[0] << "\t" << error[1] <<  endl; 

	}
	logValidationError.close();
	///////////////
	//write tp and flag file
	if(e==(Opt.getNEpochs()-1)){
		//sort AAE
		//bubble sort
		Float* temp = new Float[2];
		for(int l=(Opt.getNEpochs()-1);l>0;l--){
			for(int m=1;m<=l;m++){
				if(AAE[m-1][1]>AAE[m][1]){
					temp[0] = AAE[m][0];
					temp[1] = AAE[m][1];
					AAE[m][0] = AAE[m-1][0];
					AAE[m][1] = AAE[m-1][1];
					AAE[m-1][0] = temp[0];
					AAE[m-1][1] = temp[1];	
				}
			}
		}
		delete[] temp;
		//////////
		//update lowestAAE and write models
		ofstream bestModels;
		bestModels.open((modelDir+"/models").c_str(),ios_base::app);
		for(int i=0;i<numberOfModels;i++){
			lowestAAE[i] = AAE[i][0];
			//writing models
			if(!bestModels.is_open()){
				cerr << "CANNOT OPEN "+modelDir+"/models" << endl;
			}	
			bestModels << AAE[i][0] << endl;
			//debug
			//cout << "Adding " << AAE[i][0] << " to models" << endl;
			///////
			////////////////
		}
		bestModels.close();
		///////////////////////////////////
		//debug
		//for(int i=0;i<(Opt.getNEpochs()-1);i++){
		//	cout << "Epoch=" << AAE[i][0] << " AAE="  << AAE[i][1] << endl; 
		//}
		//for(int i=0;i<numberOfModels;i++){
		//	cout << "Best model " << i << "Epoch= " << lowestAAE[i] << endl;  
		//}
		///////
		//remove unnecessary models
		for(int i=numberOfModels;i<Opt.getNEpochs();i++){
			stringstream es;
			es << AAE[i][0];
			if(system(("rm "+(modelDir+"/MNN"+es.str())).c_str())!=0){
				cerr << "CANNOT RUN " << "rm "+(modelDir+"/MNN"+es.str()) << endl;
				exit(0);	
			}if(system(("rm "+(modelDir+"/ONN"+es.str())).c_str())!=0){
				cerr << "CANNOT RUN " << "rm "+(modelDir+"/ONN"+es.str()) << endl;
				exit(0);	
			}

		}	
		//////////////////////////
		//compute bestP
		vector<Float> modelOut(validationVector.size(),0); 
		for(int i=0;i<validationVector.size();i++){
			for(int j=0;j<numberOfModels;j++){
				//read models
				stringstream es;
				es << lowestAAE[j];
				//debug
				//cout << modelDir+"/ONN"+es.str() << endl;
				///////
				read(modelDir+"/ONN"+es.str(),modelDir+"/MNN"+es.str());
				/////////////
				dummy[0] = validation[i]->logP;
				modelOut[i]+=forward(validation[i],dummy);
			}
			bestP[i] = denormalize(modelOut[i]/(Float)numberOfModels,bound);
			//bestP[i] = modelOut[i]/(Float)numberOfModels;
			//bestT[i] = denormalize(validation[i]->target,bound);	
			bestT[i] = validation[i]->trueTarget;
				
		}	
		///////////////
		//debug
		cout << "Writing logTPVal" << endl;
		///////
		for(int i=0;i<validationVector.size();i++){
			logTPVal << bestP[i] << "\t" << bestT[i] << endl;
			//debug
			cout << bestP[i] << "\t" << bestT[i] << endl;
			///////
		}
		ofstream flagFile;
		flagFile.open((modelDir+"/flagVal").c_str());
		if(!flagFile.is_open()){
			cerr << "CANNOT OPEN FLAG FILE" << endl;
			exit(1);
		}
		flagFile << "DONE";
		//flagFile.close();
		logTPVal.close();
	}
	////////////////////////
	//free memory
	delete [] dummy;
	/////////////	
}

void UGRNN::predict(Molecule** test, const vector<Float> & bound,int numberOfModels){
	//predict targets 
	vector<Float> modelOut(testVector.size(),0); 
	vector<string> models(numberOfModels,"");
	vector<Float> p(testVector.size(),0);
	vector<Float> t(testVector.size(),0);
	Float* dummy = new Float[4];
	memset(dummy,0.0,sizeof(Float)*4);
	//read list of models from file
	ifstream modelsFile;
	modelsFile.open((modelDir+"/models").c_str());
	
	if(!modelsFile.is_open()){
		cerr << "CANNOT OPEN models FILE" << endl;
		exit(1);
	}
	for(int i=0;i<numberOfModels;i++){
		modelsFile >> models[i];	
	}	
	//debug
	for(int i=0;i<numberOfModels;i++){
		cout << models[i] << endl;
	}
	///////	
	///////////////////////////////
	for(int i=0;i<testVector.size();i++){
		for(int j=0;j<numberOfModels;j++){
			//read models
			//debug
			cout << modelDir+"/ONN"+models[j] << endl;
			///////
			read(modelDir+"/ONN"+models[j],modelDir+"/MNN"+models[j]);
			/////////////
			dummy[0] = (Float)test[i]->logP;
			modelOut[i]+=forward(test[i],dummy);
		}
		p[i] = denormalize(modelOut[i]/(Float)numberOfModels,bound);
		t[i] = test[i]->trueTarget;
	}	
	///////////////
	//debug
	cout << "Writing logTP" << endl;
	///////
	for(int i=0;i<testVector.size();i++){
		logTP << p[i] << "\t" << t[i] << endl;
		//debug
		cout << p[i] << "\t" << t[i] << endl;
		///////
	}
	ofstream flagFile;
	flagFile.open((modelDir+"/flag").c_str());
	if(!flagFile.is_open()){
		cerr << "CANNOT OPEN FLAG FILE" << endl;
		exit(1);
	}
	flagFile << "DONE";
	//flagFile.close();
	logTP.close();
		
	//free memory
	delete[] dummy;
	/////////////

}


vector<Float> UGRNN::predictTarget(Molecule** input, const vector<Float> bound, int numberOfModels){
	vector<Float> modelOut(inputVector.size(),0); 
	vector<string> models(numberOfModels,"");
	vector<Float> p(inputVector.size(),0);
	Float* dummy = new Float[1];
	memset(dummy,0.0,sizeof(Float));
	//read list of models from file
	ifstream modelsFile;
	modelsFile.open((modelDir+"/models").c_str());
	if(!modelsFile.is_open()){
		cerr << "CANNOT OPEN models FILE" << endl;
		exit(1);
	}
	for(int i=0;i<numberOfModels;i++){
		modelsFile >> models[i];	
	}	
	//debug
	//for(int i=0;i<numberOfModels;i++){
	//	cout << models[i] << endl;
	//}
	//cout << "inputVector.size()=" << inputVector.size() << endl;
	///////	
	///////////////////////////////
	
	for(int i=0;i<inputVector.size();i++){
		for(int j=0;j<numberOfModels;j++){
			Float out;			
			//read models
			//debug
			cout << modelDir+"/ONN"+models[j] << endl;
			///////
			read(modelDir+"/ONN"+models[j],modelDir+"/MNN"+models[j]);
			/////////////
			dummy[0] = (Float)input[i]->logP;
			//the following statement should be updated with a switch-case             
			if(propertiesNum == 1){			
                dummy[0] = normalizeLogP(dummy[0],bound);
            }
			out = forward(input[i],dummy);			
			modelOut[i]+=out;
			//debug
			cout << denormalize(out,bound) << endl;
			///////
		}
		p[i] = denormalize(modelOut[i]/(Float)numberOfModels,bound);
	}	
	///////////////
	//free memory
	delete[] dummy;
	/////////////	
	
	return p;		
}

void UGRNN::printDescriptors(Molecule** train,Molecule** test,vector<Float> bound,int numberOfModels){
	//global variables
	vector<string> models(numberOfModels,"");
	Float* dummy = new Float[1];
	memset(dummy,0,sizeof(Float));
	Float* descriptors = new Float[Opt.getMNN_OUTPUTS_NUMBER()];
	//////////////////
	//read list of models from file
	ifstream modelsFile;
	modelsFile.open((modelDir+"/models").c_str());
	if(!modelsFile.is_open()){
		cerr << "CANNOT OPEN models FILE" << endl;
		exit(1);
	}
	for(int i=0;i<numberOfModels;i++){
		modelsFile >> models[i];	
	}	
	//debug
	for(int i=0;i<numberOfModels;i++){
		cout << models[i] << endl;
	}
	///////	
	////////////////////////////////
	Float modelOut = 0;
	Float target = 0;
	for(int j=0;j<numberOfModels;j++){
		//create trainDescriptor and testDescriptor files
		setLogTrainDescriptors(modelDir+"/trainDescriptors"+models[j]);
		setLogTestDescriptors(modelDir+"/testDescriptors"+models[j]);					
		///////////////////////////////////////////////// 
		for(int i=0;i<testVector.size();i++){
			//read models
			//debug
			//cout << modelDir+"/ONN"+models[j] << endl;
			///////
			read(modelDir+"/ONN"+models[j],modelDir+"/MNN"+models[j]);
			/////////////
			dummy[0] = (Float)test[i]->logP;
			memset(descriptors,0,sizeof(Float)*Opt.getMNN_OUTPUTS_NUMBER());
			modelOut = forward(test[i],dummy,descriptors);
			//debug
			//for(int d=0;d<MNN->get_NO();d++){
			//	cout << descriptors[d] << " ";
			//}
			//cout << endl;
			///////
			modelOut = denormalize(modelOut,bound);
			target = test[i]->trueTarget;
			//write log file
			logTestDescriptors << testVector[i] << "\t";
			for(int d=0;d<MNN->get_NO();d++){
				logTestDescriptors << descriptors[d] << "\t";
			}
			logTestDescriptors << modelOut << "\t" << target << endl;	
			////////////////
		}
		for(int i=0;i<trainVector.size();i++){
			//read models
			//debug
			//cout << modelDir+"/ONN"+models[j] << endl;
			///////
			read(modelDir+"/ONN"+models[j],modelDir+"/MNN"+models[j]);
			/////////////
			dummy[0] = (Float)train[i]->logP;
			memset(descriptors,0,sizeof(Float)*Opt.getMNN_OUTPUTS_NUMBER());
			modelOut = forward(train[i],dummy,descriptors);
			//debug
			//for(int d=0;d<MNN->get_NO();d++){
			//	cout << descriptors[d] << " ";
			//}
			//cout << endl;
			///////
			modelOut = denormalize(modelOut,bound);
			target = train[i]->trueTarget;
			//write log file
			logTrainDescriptors << trainVector[i] << "\t";
			for(int d=0;d<MNN->get_NO();d++){
				logTrainDescriptors << descriptors[d] << "\t";
			}
			logTrainDescriptors << modelOut << "\t" << target << endl;	
			////////////////
		}
		//close logFiles
		logTrainDescriptors.close();
		logTestDescriptors.close();
		/////////////
	}
	///////////////
	/////////////////
	//free memory
	delete[] descriptors;
	delete[] dummy;
}

Float* UGRNN::compute(int node,NodeIO** nodeIO,int** children,int* num_children,int numChildren,int** bond_type){
	/*	
	for(int a=0;a<26;a++){
		cout << "chidlren for node " << a << endl;
		for(int b=0;b<num_children[a];b++){
			cout << children[a][b] << " ";
		}
		cout << endl;
	}
	*/		
	

	//cout << "node " << node << endl;
	//cout << "num children " << numChildren << endl;
	//@debug
	//cout << "computing node: " << node << endl;
	//@debug 
	//cout << "label:";
	//nodeIO[node]->printI1();
	//cout << "numChildren: " << numChildren << endl;
	if(numChildren==0){
		/*
		cout << "leaf" << endl;
                cout << "forwarding" << node << endl;
                nodeIO[node]->printI1();
                nodeIO[node]->printI2();
		*/
		//@debug
		//cout << "forwarding node: " << node << endl; 
		//nodeIO[node]->printI1();
                //nodeIO[node]->printI2();
                MNN->forward(nodeIO[node]->getI1(),nodeIO[node]->getI2());
                nodeIO[node]->setOut(MNN->out());
                return MNN->out();
	
	}else{
		//@debug
		//for(int i=0;i<numChildren;i++){
                //	cout << children[node][i] << " " << endl;
		//}
		//cout << endl;
		int b = nodeIO[node]->getOutputSize();
		int c = 0;
		for(int i=0;i<numChildren;i++){
			nodeIO[children[node][i]]->setOut(compute(children[node][i],nodeIO,children,num_children,num_children[children[node][i]],bond_type));
                	//nodeIO[children[node][i]]->printOut();
			for(int a=c;a<b;a++){ 
                        	nodeIO[node]->setI2(a,nodeIO[children[node][i]]->getOut(a-c));
                	}
			nodeIO[node]->setI2(b,bond_type[node][children[node][i]]);
			b+=nodeIO[node]->getOutputSize()+1;
			c+=nodeIO[node]->getOutputSize()+1;	
		}
		//cout << "print out for node:" << node << endl;	
		//debug routine
		//nodeIO[node]->printI2();
		/*
		cout << "forwarding " << node << endl;
                nodeIO[node]->printI1();
                nodeIO[node]->printI2();
		*/
		//cout << "forwarding node: " << node << endl;
		//nodeIO[node]->printI1();
                //nodeIO[node]->printI2();
		MNN->forward(nodeIO[node]->getI1(),nodeIO[node]->getI2());
                nodeIO[node]->setOut(MNN->out());
                return MNN->out();	
	}
		
}
void UGRNN::computeBack(int node,NodeIO** nodeIO,Float* delta,int** children,int* num_children,int numChildren){
	//@debug
	//cout << "computing back for node " << node << endl;
	//cout << "numChildren: " << numChildren << endl;
	//for(int i=0;i<numChildren;i++){
	//	cout << children[node][i] << " ";
	//}
	//cout << endl;
	///////	
	MNN->forward(nodeIO[node]->getI1(),nodeIO[node]->getI2());
	//cout << "delta=" << delta << endl; 
	//@debug
	//cout << "input delta: " << endl;
	//for(int i=0;i<(nodeIO[node]->getOutputSize());i++){
	//	cout << delta[i] << " ";
	//}
	//cout << endl;
	////////
	MNN->backward(delta);
	MNN->gradient(nodeIO[node]->getI1(),nodeIO[node]->getI2(),delta);
	if(numChildren!=0){
		Float* currentDelta = new Float[(MAX_CHILDREN)*(nodeIO[node]->getOutputSize()+1)];
		memset(currentDelta,0,sizeof(Float)*(MAX_CHILDREN)*(nodeIO[node]->getOutputSize()+1));
		//cout << "ciao" << endl;
		for(int j=0;j<numChildren;j++){
			//debug
			//cout << "delta out: " << endl; 
			///////
			for(int i=0;i<nodeIO[node]->getOutputSize();i++){
				currentDelta[i+j*(nodeIO[node]->getOutputSize()+1)] = MNN->back_out()[ATOM_NUMBER+i+j*(nodeIO[node]->getOutputSize()+1)]; 
				//@debug
				//cout << currentDelta[i+j*(nodeIO[node]->getOutputSize()+1)] << " ";
				////////
			}
			//@debug
			//cout << endl;
			////////
		}
		for(int j=0;j<numChildren;j++){
			Float* deltaOut = new Float[nodeIO[node]->getOutputSize()];	
			memset(deltaOut,0,sizeof(Float)*(nodeIO[node]->getOutputSize()));
                        //@debug
			//cout << "delta for children: " << children[node][j] << " with father: "  << node << endl;
			//////// 
			for(int i=0;i<(nodeIO[node]->getOutputSize());i++){
				deltaOut[i] = currentDelta[i+j*(nodeIO[node]->getOutputSize()+1)];
				//@debug
				//cout << deltaOut[i] << " ";
				////////
			}
                        //@debug
			//cout << endl;
			////////  
			computeBack(children[node][j],nodeIO,deltaOut,children,num_children,num_children[children[node][j]]);
			delete[] deltaOut;
		}
		delete[] currentDelta;
	}
	
}

//add outGraph to out
void UGRNN::sum(Float* out,Float* outSubGraph,int outputNumber){
	//cout << "outputNumber " << outputNumber << endl;
	for(int i=0;i<outputNumber;i++){
		//cout << "output " << i << endl; 
		//cout << outSubGraph[i] << endl;
		out[i]+=outSubGraph[i];
	}
}

Float UGRNN::computeAverageAtomsNumber(Molecule** moleculeArray,int moleculeNumber){
	Float averageAtomsNumber = 0;
	for(int i=0;i<moleculeNumber;i++){
		averageAtomsNumber += moleculeArray[i]->g->total_vert;
	}
	return (averageAtomsNumber/(Float)moleculeNumber);
}

void UGRNN::normalizeOut(Float* out,int outputNumber,int totalAtoms){
	for(int i=0;i<outputNumber;i++){
		//out[i] = out[i]/(Float)2;
		//out[i] = out[i]/(Float)outputNumber;
		out[i] = (Float)out[i]/(Float)totalAtoms;
		//out[i] = 0.01*out[i];
	}
}

void UGRNN::computeAverageOut(Float* totalOut,int outputNumber,int moleculeNumber){
	for(int i=0;i<outputNumber;i++){
		totalOut[i]=totalOut[i]/(Float)moleculeNumber;
	}
}

void UGRNN::print(Float* array,int outputNumber,string message){
	cout << message;
	for(int i=0;i<outputNumber;i++){
		cout << array[i] << " ";
	}
	cout << endl;
}

string UGRNN::unTranslateU(int a){
        if(a==0) return "C";
        else if(a==1) return "H";
        else if(a==2) return "N";
        else if(a==3) return "O";
        else if(a==4) return "S";
        else if(a==5) return "F";
        else if(a==6) return "Cl";
        else if(a==7) return "Br";
        else if(a==8) return "I";
        else if(a==9) return "In";
        else if(a==10) return "K";
        else if(a==11) return "Na";
        else if(a==12) return "Ba";
        else if(a==13) return "Sb";
        else if(a==14) return "P";
        else if(a==15) return "Be";
        else if(a==16) return "Sn";
        else if(a==17) return "Cu";
        else if(a==18) return "B";
        else if(a==19) return "Cd";
        else if(a==20) return "Ca";
        else if(a==21) return "As";
        else if(a==22) return "Co";
        else if(a==23) return "Cr";
        else if(a==24) return "Te";
        else if(a==25) return "Fe";
        else if(a==26) return "Pb";
        else if(a==27) return "Mn";
        else if(a==28) return "Hg";
        else if(a==29) return "Mo";
        else if(a==30) return "Ni";
        else if(a==31) return "Se";
        else if(a==32) return "Ti";
        else if(a==33) return "Zn";
        else if(a==34) return "Si";
        else if(a==35) return "?";
        else if(a==36) return "Mg";
        else if(a==37) return "V";
        else if(a==38) return "Li";
        else if(a==39) return "Al";
        else if(a==40) return "Zr";
        else if(a==41) return "Bi";
        else if(a==42) return "Pd";
	else if(a==43) return "Pt";
        else if(a==44) return "Ru";
        else if(a==45) return "Rh";
        else if(a==46) return "Ga";
        else if(a==47) return "Ge";
        else if(a==48) return "Ag";
        else if(a==49) return "Tb";
        else if(a==50) return "Ir";
        else if(a==51) return "W";
        else if(a==52) return "Cs";
        else if(a==53) return "Re";
        else if(a==54) return "Pr";
        else if(a==55) return "Nd";
        else if(a==56) return "Gd";
        else if(a==57) return "Yb";
        else if(a==58) return "Er";
        else if(a==59) return "U";
        else if(a==60) return "Tl";
        else if(a==61) return "Au";
        else if(a==62) return "Ac";
        else if(a==63) return "Ho";
        else if(a==64) return "Os";
        else if(a==65) return "Sm";
        else if(a==66) return "Nb";
        else if(a==67) return "R3";
        else if(a==68) return "R4";
        else if(a==69) return "R5";
        else if(a==70) return "R6";
        else if(a==71) return "R7";
        else if(a==72) return "R8";
        else if(a==73) return "R14";
        else if(a==74) return "R24";
        else return "None";

}
vector<string> UGRNN::getAllAtoms(Molecule** moleculeList,int moleculeNumber){
        vector<string> atomSymbol;
        int counter;
        string symbol="";
        for(int i=0;i<moleculeNumber;i++){
                for(int j=0;j<moleculeList[i]->atoms;j++){
                        counter = 0;
                        symbol=UGRNN::unTranslateU(moleculeList[i]->g->vert[j]->symbol);
                        //cout << symbol;
                        for(int k=0;k<atomSymbol.size();k++){
                                if(symbol.compare(atomSymbol[k])==0){
                                        counter +=1;
                                }
                        }
                        if (counter!=1){
                                atomSymbol.push_back(symbol);
                        }
                }

            
        }
        return atomSymbol;
}
int UGRNN::getMaxNumberOfChildren(Molecule** moleculeList,int moleculeNumber){
	int max = 0;
	for(int i=0;i<moleculeNumber;i++){
		moleculeList[i]->g->orderGraph();
		moleculeList[i]->g->loadNeighbours();		
		for(int j=0;j<moleculeList[i]->g->total_vert;j++){
			for(int k=0;k<moleculeList[i]->g->total_vert;k++){
				if(max<moleculeList[i]->g->num_children[j][k]){
					max = moleculeList[i]->g->num_children[j][k];
				}

			}	
		}
		moleculeList[i]->g->unloadGraph();	
	}	
	return max;
}

vector<int> UGRNN::shuffle(int moleculeNumber,int seed){
	vector<int> out(moleculeNumber,-1);
	srand(seed);
	
	int tempPos;
	for(int i=0;i<moleculeNumber;i++){
		//cout << "shuffle molecule " << i << endl;
		tempPos = (int)(moleculeNumber*(rand()/(RAND_MAX+1.0)));
		while(out[tempPos]!=-1){
			tempPos = (int)(moleculeNumber*(rand()/(RAND_MAX+1.0)));
		}
		out[tempPos] = i;
		//cout << i << " " << tempPos << endl;
	}
	return out;
}

void UGRNN::initWeights(int seed){
	MNN->initWeights(seed);
	ONN->initWeights(seed);
}

//bound[0]==min, bound[1]=max, bound[2]=myMin, bound[3]=myMax
Float UGRNN::denormalize(Float value,const vector<Float> & bound){
	Float out = 0;
	out = (Float)(value-bound[2])*(Float)(bound[1]-bound[0])/(Float)(bound[3]-bound[2])+(Float)bound[0];	
	//(Float)(value-myMin)*(Float)(max-min)/(Float)(myMax-myMin)+(Float)min;	
	return out;
}

//bound[0]==min, bound[1]=max, bound[2]=myMin, bound[3]=myMax, bound[4]=minLogP, bound[5]=maxLogP
Float UGRNN::normalizeLogP(Float value,const vector<Float> & bound){
	//(Float)myMin+(Float)myRange*(molecule[i]->logP / (Float)(maxLogP-minLogP));	
	return (Float)bound[2]+(Float)(bound[3]-bound[2])*(value / (Float)(bound[5]-bound[4]));
}

Float UGRNN::computeCorrelation(const vector<Float> & t,const vector<Float> & p,int N){
	Float r = 0;
	//average t
	Float aT = 0;
	//average p
	Float aP = 0;
	//st
	Float st = 0;
	//sp
	Float sp = 0;
	//calculating average t and p
	for(int i=0;i<N;i++){
		aT += t[i]; 
		aP += p[i];
	}
	aT = aT/(Float)N;
	aP = aP/(Float)N;
	///////////////////////
	//calculating st and sp
	for(int i=0;i<N;i++){
		st += (t[i]-aT)*(t[i]-aT);
		sp += (p[i]-aP)*(p[i]-aP);
	}
	///////////////////////
	//calculating r
	for(int i=0;i<N;i++){
		r += t[i]*p[i];
	}
	r = (r-(Float)N*aT*aP)/sqrt(st*sp);
			
	return (r*r);
}

vector<Float> UGRNN::computeError(const vector<Float> & t,const vector<Float> & p,int N){
	vector<Float> out(2,0);
	//compute RMSE
	for(int i=0;i<N;i++){
		out[0] = out[0] + (t[i]-p[i])*(t[i]-p[i]);
	}
	out[0] = (Float)sqrt(out[0]/(Float)N);
	//////////////
	//compute AAE
	for(int i=0;i<N;i++){
		out[1] = out[1] + fabs((t[i]-p[i]));
	}
	out[1] = (Float)(out[1]/(Float)N);
	//////////////
	return out;
}

void UGRNN::writeTotalAtoms(string dir){
	ofstream out;
	//debug
	cout << dir << endl;
	///////
	out.open(dir.c_str());
	if(out.is_open()){
		for(int i=0;i<totalAtoms.size()-1;i++){
			out << totalAtoms[i] << endl;
		}
		out << totalAtoms[totalAtoms.size()-1] << endl;
		out.close();
	}else{
		cerr << "CANNOT OPEN totalAtoms FILE" << endl;
		exit(1);
	}
}

void UGRNN::readTotalAtoms(string dir){
	ifstream in;
	in.open(dir.c_str());
	if(in.is_open()){
		string atom;	
		in >> atom;
		while(in){
			//debug
			//cout << atom << endl;
			///////
			totalAtoms.push_back(atom);
			in >> atom;
		}	
		in.close();
	}else{
		cerr << "CANNOT OPEN totalAtoms FILE" << endl;
	}
}

void UGRNN::compareTotalAtoms(vector<string> totalAtoms,vector<string> totalAtomsTest){
	vector<int> temp1;
	vector<int> temp2;
	
	//debug
	//cout << "totalAtoms_size " << totalAtoms.size() << endl;
	//for(int i=0;i<totalAtoms.size();i++){
	//	cout << totalAtoms[i] << " ";
	//}
	//cout << endl;
	//cout << "totalAtoms_size " << totalAtomsTest.size() << endl;
	//for(int i=0;i<totalAtomsTest.size();i++){
	//	cout << totalAtomsTest[i] << " ";
	//}
	//cout << endl;
	///////
	//Convert element types into numbers
	for(int i=0;i<totalAtoms.size();i++){
		temp1.push_back(Molecule::translateU(totalAtoms[i]));		
	}
	for(int i=0;i<totalAtomsTest.size();i++){
		temp2.push_back(Molecule::translateU(totalAtomsTest[i]));
	} 	
	////////////////////////////////////
	//sort vectors
	sort(temp1.begin(),temp1.end());
	sort(temp2.begin(),temp2.end());
	//////////////
	
	//compare vectors
	int check = 0;
	for(int i=0;i<temp2.size();i++){		
		for(int j=0;j<temp1.size();j++){
			if(temp1.at(j) == temp2.at(i)){
				check = 1;	
			}
		}
		if(check == 0){
				cerr <<"TEST SET CONTAINS AN ATOM WHICH IS NOT INSIDE THE TRAINING SET" << endl;
				exit(0);		
	   }else{
			check = 0;		
	   }
	}	
	//exit(0);	
	/////////////////
}

void UGRNN::writeMAX_CHILDREN(string dir){
	ofstream out;
	out.open(dir.c_str());
	if(out.is_open()){
		out << this->MAX_CHILDREN;
	}else{
		cerr << "CANNOT OPEN MAX_CHILDREN FILE FOR WRITING" << endl;
	}
	out.close();	
}

void UGRNN::readMAX_CHILDREN(string dir){
	ifstream in;
	in.open(dir.c_str());
	if(in.is_open()){
		in >> this->MAX_CHILDREN;
	}else{
		//debug
		cerr << dir << "endl";
		///////
		cerr << "CANNOT OPEN MAX_CHILDREN FILE FOR READING" << endl;
		exit(0);
	}
	in.close();
}
