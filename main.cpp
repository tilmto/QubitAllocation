#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#define infinity 1e4

using namespace std;

///////////////  Beginning of Hardware Description  ///////////////////

class Hardware
{
private:
    int phyQubitNum;

    int edgeNum;

    bool isUniDirection;

    vector<vector<bool>> archMatrix;

    vector<vector<int>> distMatrix;

public:
    vector<vector<int>> routeMatrix;

    vector<int> outdeg;

    Hardware(string hwname,bool isUniDirection);

    int GetQNum();

    int GetENum();

    void VerifyArchMatrix();

    void PrintArchMatrix();

    void Floyd();

    void VerifyRouteMatrix();

    void PrintRouteMatrix();

    void PrintPath(int i,int j);
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

    phyQubitNum=0;
    edgeNum=0;

    while(!is.eof())
    {
        is>>adjIndex;
        if(adjIndex == -1)
            phyQubitNum++;
    }
    phyQubitNum--;

    for(i=0; i<phyQubitNum; i++)
    {
        archMatrix.push_back(vector<bool>(phyQubitNum,false));
        distMatrix.push_back(vector<int>(phyQubitNum));
        routeMatrix.push_back(vector<int>(phyQubitNum));
        outdeg.push_back(0);
    }

    i=0;
    is.clear();
    is.seekg(0,ios::beg);

    while(i<phyQubitNum && !is.eof())
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

    cout << "Physical qubits number: " << phyQubitNum << endl;
    cout << "Edge number: " << edgeNum << endl;
}


void Hardware::VerifyArchMatrix()
{
    for(int i=0; i<phyQubitNum; i++)
        for(int j=i; j<phyQubitNum; j++)
            if(archMatrix[i][j]!=archMatrix[j][i])
            {
                cout << "Wrong Hardware Architecture." << endl;
                exit(1);
            }
}


void Hardware::PrintArchMatrix()
{
    cout << "Architecture Matrix:" << endl;
    for(int i=0; i<phyQubitNum; i++)
        for(int j=0; j<phyQubitNum; j++)
        {
            cout << archMatrix[i][j] << " ";
            if(j==phyQubitNum-1)
                cout << endl;
        }
}


void Hardware::Floyd()
{
    int i,j,k;
    for(i=0; i<phyQubitNum; i++)
        for(j=0; j<phyQubitNum; j++)
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

    for(k=0; k<phyQubitNum; k++)
        for(i=0; i<phyQubitNum; i++)
            for(j=0; j<phyQubitNum; j++)
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
    for(int i=0; i<phyQubitNum; i++)
        for(int j=0; j<phyQubitNum; j++)
            if(routeMatrix[i][j]==-1)
            {
                cout << "Not fully connected architecture." << endl;
                exit(1);
            }
}


void Hardware::PrintRouteMatrix()
{
    cout << "Route Matrix:" << endl;
    for(int i=0; i<phyQubitNum; i++)
        for(int j=0; j<phyQubitNum; j++)
        {
            cout << routeMatrix[i][j] << " ";
            if(j==phyQubitNum-1)
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
    return phyQubitNum;
}


int Hardware::GetENum()
{
    return edgeNum;
}

///////////////  Ending of Hardware Description  ///////////////////


void SeqGenerate(vector<vector<int>> &seq,int qubitNum,int seqLen);

void PrintSeq(vector<vector<int>> seq);

void PrintMap(vector<int> mapArray);

void InitMap(vector<int> &mapArray,Hardware hdware,vector<vector<int>> seq);

int AllocBySwap(vector<int> &mapArray,Hardware hdware,vector<vector<int>> seq);

///////////////  Main Function  ///////////////////

int main()
{
    int cost;
    Hardware hdware("ibmqx4");
    vector<int> mapArray(hdware.GetQNum());
    vector<vector<int>> seq;

    SeqGenerate(seq,hdware.GetQNum(),100);
    PrintSeq(seq);
    InitMap(mapArray,hdware,seq);
    cost=AllocBySwap(mapArray,hdware,seq);

    cout << "The total cost is: " << cost << endl;

    return 0;
}

///////////////  Main Function  ///////////////////

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


void InitMap(vector<int> &mapArray,Hardware hdware,vector<vector<int>> seq)
{
    int i;
    unsigned int j;
    int qubitNum=hdware.GetQNum();
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
            if(hdware.outdeg[i]>hdware.outdeg[sortOutDeg[j]])
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
    PrintMap(mapArray);
}


void PrintMap(vector<int> mapArray)
{
    unsigned int i;
    cout << "Physical qubits: ";
    for(i=0; i<mapArray.size(); i++)
        cout << i << " ";
    cout << endl;
    cout << "Pseudo   qubits: ";
    for(i=0; i<mapArray.size(); i++)
        cout << mapArray[i] << " ";
    cout << endl;
}


int AllocBySwap(vector<int> &mapArray,Hardware hdware,vector<vector<int>> seq)
{
    unsigned int i;
    int j,temp,current,next,dest,cost=0;
    int qubitNum=hdware.GetQNum();

    for(i=0; i<seq.size(); i++)
    {
        for(j=0; j<qubitNum; j++)
        {
            if(mapArray[j]==seq[i][1])
                current=j;

            if(mapArray[j]==seq[i][0])
                dest=j;
        }

        next=hdware.routeMatrix[current][dest];

        while(next!=dest)
        {
            temp=mapArray[current];
            mapArray[current]=mapArray[next];
            mapArray[next]=temp;
            cost=cost+7;
            current=next;
            next=hdware.routeMatrix[current][dest];
        }

        cost++;

        cout << "After the " << i << "th operation:" << endl;
        PrintMap(mapArray);
        cout << endl;
    }

    return cost;
}

















