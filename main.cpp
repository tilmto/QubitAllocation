#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#define infinity 1e4

using namespace std;

///////////////  Beginning of Original Hardware Description  ///////////////////

class Hardware
{
protected:
    int qubitNum;

    int edgeNum;

    bool isUniDirection;

    vector<vector<bool>> archMatrix;

    vector<vector<int>> distMatrix;

    vector<vector<int>> routeMatrix;

    vector<int> outdeg;

    vector<int> mapArray;

public:
    Hardware(string hwname,bool isUniDirection);

    int GetQNum();

    int GetENum();

    void VerifyArchMatrix();

    void PrintArchMatrix();

    void Floyd();

    void VerifyRouteMatrix();

    void PrintRouteMatrix();

    void PrintPath(int i,int j);

    void PrintMap();

    void InitMap(vector<vector<int>> seq);

    int Alloc(vector<vector<int>> seq);

};


Hardware::Hardware(string hwname,bool isUniDirection=false)
{
    int adjIndex,i;

    this->isUniDirection=isUniDirection;

    ifstream is(hwname,ios::in);
    if(!is)
    {
        cout << "Cannot Open Hardware File." << endl;
        exit(1);
    }

    qubitNum=0;
    edgeNum=0;

    while(!is.eof())
    {
        is>>adjIndex;
        if(adjIndex == -1)
            qubitNum++;
    }
    qubitNum--;

    for(i=0; i<qubitNum; i++)
    {
        archMatrix.push_back(vector<bool>(qubitNum,false));
        distMatrix.push_back(vector<int>(qubitNum));
        routeMatrix.push_back(vector<int>(qubitNum));
        outdeg.push_back(0);
    }

    mapArray.resize(qubitNum);

    i=0;
    is.clear();
    is.seekg(0,ios::beg);

    while(i<qubitNum && !is.eof())
    {
        is>>adjIndex;
        if(adjIndex==-1)
        {
            archMatrix[i][i]=true;
            i++;
        }

        else
        {
            archMatrix[i][adjIndex]=true;
            outdeg[i]++;
            edgeNum++;
        }
    }

    is.close();

    edgeNum=edgeNum/2;

    VerifyArchMatrix();

    PrintArchMatrix();

    Floyd();

    cout << "Physical qubits number: " << qubitNum << endl;
    cout << "Edge number: " << edgeNum << endl;
}


void Hardware::VerifyArchMatrix()
{
    for(int i=0; i<qubitNum; i++)
        for(int j=i; j<qubitNum; j++)
            if(archMatrix[i][j]!=archMatrix[j][i])
            {
                cout << "Wrong Hardware Architecture." << endl;
                exit(1);
            }
}


void Hardware::PrintArchMatrix()
{
    cout << "Architecture Matrix:" << endl;
    for(int i=0; i<qubitNum; i++)
        for(int j=0; j<qubitNum; j++)
        {
            cout << archMatrix[i][j] << " ";
            if(j==qubitNum-1)
                cout << endl;
        }
}


void Hardware::Floyd()
{
    int i,j,k;
    for(i=0; i<qubitNum; i++)
        for(j=0; j<qubitNum; j++)
        {
            if(!archMatrix[i][j])
            {
                distMatrix[i][j]=infinity;
                routeMatrix[i][j]=-1;
            }

            else
            {
                distMatrix[i][j]=1;
                routeMatrix[i][j]=j;
            }
        }

    for(k=0; k<qubitNum; k++)
        for(i=0; i<qubitNum; i++)
            for(j=0; j<qubitNum; j++)
                if(distMatrix[i][j]>distMatrix[i][k]+distMatrix[k][j])
                {
                    distMatrix[i][j]=distMatrix[i][k]+distMatrix[k][j];
                    routeMatrix[i][j]=routeMatrix[i][k];
                }

    VerifyRouteMatrix();

    PrintRouteMatrix();
}



void Hardware::VerifyRouteMatrix()
{
    for(int i=0; i<qubitNum; i++)
        for(int j=0; j<qubitNum; j++)
            if(routeMatrix[i][j]==-1)
            {
                cout << "Not fully connected architecture." << endl;
                exit(1);
            }
}


void Hardware::PrintRouteMatrix()
{
    cout << "Route Matrix:" << endl;
    for(int i=0; i<qubitNum; i++)
        for(int j=0; j<qubitNum; j++)
        {
            cout << routeMatrix[i][j] << " ";
            if(j==qubitNum-1)
                cout << endl;
        }
}


void Hardware::PrintPath(int i,int j)
{
    int next=routeMatrix[i][j];
    if(next==-1)
        cout << "No Path between " << i << " and "<< j << endl;
    else
    {
        cout << "Path from " << i << " to " << j << ": " << i << " ";
        while(next!=j)
        {
            cout << next << " ";
            next=routeMatrix[next][j];
        }
        cout << j << endl;
    }
}

int Hardware::GetQNum()
{
    return qubitNum;
}


int Hardware::GetENum()
{
    return edgeNum;
}


void Hardware::InitMap(vector<vector<int>> seq)
{
    int i;
    unsigned int j;
    vector<int> freq(qubitNum,0);
    vector<int> sortFreq(1,0);
    vector<int> sortOutDeg(1,0);

    for(j=0; j<seq.size(); j++)
        freq[seq[j][0]]++;

    for(i=1; i<qubitNum; i++)
        for(j=0; j<sortFreq.size(); j++)
        {
            if(freq[i]>freq[sortFreq[j]])
            {
                sortFreq.insert(sortFreq.begin()+j,i);
                break;
            }

            if(j==sortFreq.size()-1)
            {
                sortFreq.push_back(i);
                break;
            }
        }

    for(i=1; i<qubitNum; i++)
        for(j=0; j<sortOutDeg.size(); j++)
        {
            if(outdeg[i]>outdeg[sortOutDeg[j]])
            {
                sortOutDeg.insert(sortOutDeg.begin()+j,i);
                break;
            }

            if(j==sortOutDeg.size()-1)
            {
                sortOutDeg.push_back(i);
                break;
            }
        }

    for(i=0; i<qubitNum; i++)
        mapArray[sortOutDeg[i]]=sortFreq[i];

    cout << "Initial Mapping:" << endl;
    PrintMap();
}


void Hardware::PrintMap()
{
    int i;
    cout << "Physical qubits: ";
    for(i=0; i<qubitNum; i++)
        cout << i << " ";
    cout << endl;
    cout << "Pseudo   qubits: ";
    for(i=0; i<qubitNum; i++)
        cout << mapArray[i] << " ";
    cout << endl;
}


int Hardware::Alloc(vector<vector<int>> seq)
{
    unsigned int i;
    int j,temp,current,next,dest,cost=0;

    for(i=0; i<seq.size(); i++)
    {
        for(j=0; j<qubitNum; j++)
        {
            if(mapArray[j]==seq[i][1])
                current=j;

            if(mapArray[j]==seq[i][0])
                dest=j;
        }

        next=routeMatrix[current][dest];

        while(next!=dest)
        {
            temp=mapArray[current];
            mapArray[current]=mapArray[next];
            mapArray[next]=temp;
            cost=cost+7;
            current=next;
            next=routeMatrix[current][dest];
        }

        cost++;

        cout << "After the " << i+1 << "th operation:" << endl;
        PrintMap();
        cout << endl;
    }

    return cost;
}

///////////////  Ending of Original Hardware Description  ///////////////////


///////////////  Beginning of New Hardware Description  /////////////////////

class HardwareNew:public Hardware
{
protected:
    vector<int> fixedQubit;

public:
    HardwareNew(string hwname,bool isUniDirection);

    void InitMap(vector<vector<int>> seq);

    int Alloc(vector<vector<int>> seq);
};

HardwareNew::HardwareNew(string hwname,bool isUniDirection=false):Hardware(hwname,isUniDirection){}

void HardwareNew::InitMap(vector<vector<int>> seq)
{
    int i;
    unsigned int j;
    vector<int> freq(qubitNum,0);
    vector<int> sortFreq(1,0);
    vector<int> sortOutDeg(1,0);

    for(j=0; j<seq.size(); j++)
    {
        freq[seq[j][0]]++;
        freq[seq[j][1]]++;
    }

    for(i=1; i<qubitNum; i++)
        for(j=0; j<sortFreq.size(); j++)
        {
            if(freq[i]>freq[sortFreq[j]])
            {
                sortFreq.insert(sortFreq.begin()+j,i);
                break;
            }

            if(j==sortFreq.size()-1)
            {
                sortFreq.push_back(i);
                break;
            }
        }

    for(i=1; i<qubitNum; i++)
        for(j=0; j<sortOutDeg.size(); j++)
        {
            if(outdeg[i]>outdeg[sortOutDeg[j]])
            {
                sortOutDeg.insert(sortOutDeg.begin()+j,i);
                break;
            }

            if(j==sortOutDeg.size()-1)
            {
                sortOutDeg.push_back(i);
                break;
            }
        }

    for(i=0; i<qubitNum; i++)
        mapArray[sortOutDeg[i]]=sortFreq[i];

    fixedQubit.push_back(6);
    fixedQubit.push_back(8);

    cout << "Initial Mapping:" << endl;
    PrintMap();
}


int HardwareNew::Alloc(vector<vector<int>> seq)
{
    unsigned int i;
    int j,temp,current,next,dest,cost=0;

    for(i=0; i<seq.size(); i++)
    {
        for(j=0; j<qubitNum; j++)
        {
            if(mapArray[j]==seq[i][1])
                current=j;

            if(mapArray[j]==seq[i][0])
                dest=j;
        }

        if(outdeg[current]>outdeg[dest])
        {
            temp=current;
            current=dest;
            dest=temp;
        }

        next=routeMatrix[current][dest];

        if(next==dest)
            cost++;
        else
        {
            while(routeMatrix[next][dest]!=dest)
            {
                if(next!=6 && next!=8)
                {
                    temp=mapArray[current];
                    mapArray[current]=mapArray[next];
                    mapArray[next]=temp;
                    cost=cost+7;
                    current=next;
                    next=routeMatrix[current][dest];
                }

                else
                {
                    for(j=0; j<qubitNum; j++)
                        if(routeMatrix[current][j]==j && routeMatrix[j][dest]!=current && j!=6 && j!=8)
                        {
                            next=j;
                            break;
                        }
                }
            }

            cost=cost+4;
        }

        cout << "After the " << i+1 << "th operation:" << endl;
        PrintMap();
        cout << endl;
    }

    return cost;
}


///////////////  Ending of New Hardware Description  ///////////////////////


void SeqGenerate(vector<vector<int>> &seq,int qubitNum,int seqLen);

void PrintSeq(vector<vector<int>> seq);


int main()
{
    int costOld,costNew;
    Hardware archOld("ibmqx4");
    HardwareNew archNew("ibmqx4");

    vector<vector<int>> seq;
    SeqGenerate(seq,archOld.GetQNum(),500);
    PrintSeq(seq);

    archOld.InitMap(seq);
    costOld=archOld.Alloc(seq);

    archNew.InitMap(seq);
    costNew=archNew.Alloc(seq);

    cout << "The total cost of old hardware is: " << costOld << endl;
    cout << "The total cost of new hardware is:" << costNew << endl;
    cout << "new cost/old cost = " << (double)costNew/costOld << endl;

    return 0;
}


void SeqGenerate(vector<vector<int>> &seq,int qubitNum,int seqLen)
{
    int cqubit,squbit;
    int i=0;
    srand((int)time(0));

    while(i<seqLen)
    {
        cqubit=rand()%qubitNum;
        squbit=rand()%qubitNum;
        if(cqubit!=squbit)
        {
            seq.push_back(vector<int>(2));
            seq[i][0]=cqubit;
            seq[i][1]=squbit;
            i++;
        }
    }
}


void PrintSeq(vector<vector<int>> seq)
{
    cout << "Dependency Sequence:"<<endl;

    for(unsigned int i=0; i<seq.size(); i++)
        cout << "( " << seq[i][0] << " , " << seq[i][1] << " )" <<endl;
}









