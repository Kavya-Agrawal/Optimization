#include<bits/stdc++.h>
using namespace std;
#define int long long

typedef pair<double, double> POINT;
typedef pair<POINT, POINT> EDGE;
#define epb edgeList.push_back
#define vpb vertices.push_back

// Declarations

const int MAX_SIZE = 1e4+1;
vector<vector<int>> grid(MAX_SIZE, vector<int>(MAX_SIZE, 0));
vector<vector<int>> psum(MAX_SIZE, vector<int>(MAX_SIZE, 0));
vector<int> starts1, ends1, starts2, ends2, maxes1, maxes2;
vector<vector<int>> dp0;
vector<vector<int>> dp1;
vector<vector<int>> dp2;

struct edgesStruct
{
    vector<EDGE> edgeList;
    vector<EDGE> upEdges;
    vector<EDGE> downEdges;
    vector<POINT> vertices;
};

void getGrid(const string& filename);
void getGridTranspose(const string& filename);
int allPositivesSum();
void makepsum();
void kadane(const vector<int>& arr, int& start, int& end, int& max_sum);
pair<int,int> getOptimalValue(int rows);
vector<int> getSubarrays(int rows, int used);
vector<int> getSubarrays(int rows, int used) ;
edgesStruct getEdges(vector<int> selected, int len);
void extend_edges(edgesStruct & es );

// End Declarations

signed main(){
    // ios::sync_with_stdio(0); cin.tie(0);
    
    int testcases;
    cout << "Enter the number of testcases : ";
    cin >> testcases;
    
    for(int tc = 1; tc <= testcases; tc++){
        cout<<endl<<"Case: "<<tc<<endl;
        string filename = "input" + string(tc < 10 ? "0" : "" ) + (to_string(tc)) + ".txt";
        getGrid(filename);
        makepsum();

        int totPos = allPositivesSum();
        int bestWidth = -1;
        int ans = 0;
        int used;
        vector<int> widths = {1, 40, 50};
        cout<< "Normal grid: " <<endl;
        for(int width : widths){
            starts1.clear();
            ends1.clear();
            starts2.clear();
            ends2.clear();
            maxes1.clear();
            maxes2.clear();
            dp0.clear();
            dp1.clear();
            dp2.clear();
            pair<int, int> temp = getOptimalValue(width);
            cout << width << " : " << (temp.first * 100.00)/totPos << "\n";
            if(temp.first > ans){
                ans = temp.first;
                used = temp.second;
                bestWidth = width;
            }
        }
        bool transpose = false;
        getGridTranspose(filename);

        cout<<endl<< "Transpose grid: "<<endl;
        for(int width : widths){
            starts1.clear();
            ends1.clear();
            starts2.clear();
            ends2.clear();
            maxes1.clear();
            maxes2.clear();
            dp0.clear();
            dp1.clear();
            dp2.clear();
            pair<int, int> temp = getOptimalValue(width);
            cout << width << " : " << (temp.first * 100.00)/totPos << "\n";
            if(temp.first > ans){
                ans = temp.first;
                used = temp.second;
                bestWidth = width;
                transpose = true;
            }
        }

        if(!transpose) getGrid(filename);

        starts1.clear();
        ends1.clear();
        starts2.clear();
        ends2.clear();
        maxes1.clear();
        maxes2.clear();
        dp0.clear();
        dp1.clear();
        dp2.clear();
        pair<int, int> temp = getOptimalValue(bestWidth);
        ans = temp.first;
        used = temp.second;

        vector<int> selected = getSubarrays(bestWidth,used);

        int veri=0;
        for(int i = 0; i < selected.size(); ++i){
            if(selected[i] == 1) veri += maxes1[i];
            else if(selected[i] == 2){
                veri += maxes1[i];
                veri += maxes2[i];
            }
        };


        cout << endl<< "Bestwidth: " << bestWidth<< ", Transpose = "<<transpose << "\n";
        cout << ans << endl;
        cout << veri << endl;

        cout<<  (ans * 100.0) / totPos<<endl;
        cout<<  (veri * 100.0) / totPos<<endl;
        
        cout << totPos << endl;
        
        edgesStruct edge_str = getEdges(selected, bestWidth);
        extend_edges(edge_str);

        int numEdges = edge_str.edgeList.size();
        
        cout << "Edges used : " << numEdges << "\n";
        
        ofstream edgeOut("57_optimization_output0" + (to_string(tc)) + ".txt");
        edgeOut << ans << "\n";
        edgeOut << numEdges << ", " << numEdges << "\n";
        for (int i = 0; i < numEdges; i++)
        {
            if(edge_str.edgeList[i].first.first > 0){
                edge_str.edgeList[i].first.first++;
            }
            if(edge_str.edgeList[i].second.first > 0){
                edge_str.edgeList[i].second.first++;
            }
            if(transpose){
                swap(edge_str.edgeList[i].first.first , edge_str.edgeList[i].first.second);
                swap(edge_str.edgeList[i].second.first , edge_str.edgeList[i].second.second);
            }
            edgeOut << "(" << edge_str.edgeList[i].first.first<<", "<<edge_str.edgeList[i].first.second<<"), ("<<
            edge_str.edgeList[i].second.first<<", "<<edge_str.edgeList[i].second.second<<")\n";

        }
    }

    return 0;
}

void getGrid(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error opening file!" << endl;
        return;
    }
    for(int i = 0; i < MAX_SIZE; ++i){
        for(int j = 0; j < MAX_SIZE; ++j){
            grid[i][j] = 0;
        }
    }
    int N;
    file >> N;
    for (int x1 = 0; x1 < N; x1++) {
        //crystals
        int x, y, value;
        file >> x >> y >> value;
        grid[x][y] += value;
    }
    int M;
    file >> M;
    for (int x1 = 0; x1 < M; x1++) {
        //void mines
        int x, y, value;
        file >> x >> y >> value;
        grid[x][y] -= value;
    }
}

void getGridTranspose(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error opening file!" << endl;
        return;
    }
    for(int i = 0; i < MAX_SIZE; ++i){
        for(int j = 0; j < MAX_SIZE; ++j){
            grid[i][j] = 0;
        }
    }
    int N;
    file >> N;
    for (int x1 = 0; x1 < N; x1++) {
        int x, y, value;
        file >> x >> y >> value;
        grid[y][x] += value;
    }
    int M;
    file >> M;
    for (int x1 = 0; x1 < M; x1++) {
        int x, y, value;
        file >> x >> y >> value;
        grid[y][x] -= value;
    }
}

int allPositivesSum() {
    int ans = 0;
    for (int i = 0; i < MAX_SIZE; ++i) {
        for (int j = 0; j < MAX_SIZE; ++j) {
            if (grid[i][j] > 0) {
                ans += grid[i][j];
            }
        }
    }
    return ans;
}

void makepsum() {
    for (int i = 0; i < MAX_SIZE; i++) {
        for (int j = 0; j < MAX_SIZE; j++) {
            psum[i][j] = grid[i][j];

            if (i > 0)
                psum[i][j] += psum[i - 1][j];
            if (j > 0)
                psum[i][j] += psum[i][j - 1];
            if (i > 0 && j > 0)
                psum[i][j] -= psum[i - 1][j - 1];
        }
    }
}

int getSum(int x1, int y1, int x2, int y2 ){
    return psum[x2][y2] - psum[x1-1][y2] - psum[x2][y1-1] + psum[x1-1][y1-1];
}

void kadane(const vector<int>& arr, int& start, int& end, int& max_sum) {
    int current_max = 0;
    max_sum = 0;
    start = 0;
    end = -1;
    int temp_start = 0;

    for (int i = 0; i < arr.size(); ++i) {
        current_max += arr[i];
        
        // Update max_sum and indices if a new maximum is found
        if (current_max > max_sum) {
            max_sum = current_max;
            start = temp_start;
            end = i;
        }

        // Reset current_max if it becomes negative
        if (current_max < 0) {
            current_max = 0;
            temp_start = i + 1;
        }
    }

    // Handle case when all elements are negative
    if (end == -1) {
        max_sum = -1e8;
        start = 0;
        end = 0;
    }
}

pair<int,int> getOptimalValue(int rows) {
    int divisions = MAX_SIZE / rows;
    dp1.assign(divisions, vector<int>(1001, 0));
    dp2.assign(divisions, vector<int>(1001, 0));
    dp0.assign(divisions, vector<int>(1001, 0));

    vector<vector<int>> subgrid(rows + 1, vector<int>(MAX_SIZE, 0));
    for (int r = 0; r <= rows; ++r) {
        for (int c = 0; c < MAX_SIZE; ++c) {
            subgrid[r][c] = grid[r][c];
        }
    }

    for (int r = 1; r <= rows; ++r) {
        for (int c = 0; c < MAX_SIZE; ++c) {
            subgrid[r][c] += subgrid[r-1][c];
        }
    }

    int start, end, max_sum;
    kadane(subgrid.back(), start, end, max_sum);
    starts1.push_back(start);
    ends1.push_back(end);
    maxes1.push_back(max_sum);

    for (int c = start; c <= end; ++c) {
        subgrid.back()[c] = INT_MIN;
    }

    kadane(subgrid.back(), start, end, max_sum);
    starts2.push_back(start);
    ends2.push_back(end);
    maxes2.push_back(max_sum);
    
    for (int i = rows + 1; i < MAX_SIZE; i += rows) {
        vector<vector<int>> subgrid(rows, vector<int>(MAX_SIZE, 0));
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < MAX_SIZE; ++c) {
                subgrid[r][c] = grid[i + r][c];
            }
        }

        for (int r = 1; r < rows; ++r) {
            for (int c = 0; c < MAX_SIZE; ++c) {
                subgrid[r][c] += subgrid[r-1][c];
            }
        }

        int start, end, max_sum;
        kadane(subgrid.back(), start, end, max_sum);
        starts1.push_back(start);
        ends1.push_back(end);
        maxes1.push_back(max_sum);

        for (int c = start; c <= end; ++c) {
            subgrid.back()[c] = INT_MIN;
        }

        kadane(subgrid.back(), start, end, max_sum);
        starts2.push_back(start);
        ends2.push_back(end);
        maxes2.push_back(max_sum);
    }

    for(int i = 0; i < divisions; i++){
        dp1[i][0] = dp1[i][1] = dp1[i][2] = dp1[i][3] = -1e8;
        dp2[i][0] = dp2[i][1] = dp2[i][2] = dp2[i][3] = dp2[i][4] = dp2[i][5] = dp2[i][6] = dp2[i][7] = -1e8;
    }
    dp1[0][4] = maxes1[0];
    dp2[0][8] = maxes1[0] + maxes2[0];

    for (int i = 1; i < divisions; ++i) {
        for (int j = 1; j <= 1000; ++j) {

            // dp0
            dp0[i][j] = max({dp0[i][j], dp0[i-1][j], dp1[i-1][j], dp2[i-1][j]});


            // dp1
            // From dp0[i - 1]
            if (j >= 6) {
                dp1[i][j] = max(dp1[i][j], dp0[i-1][j-6] + maxes1[i]);
            }

            // From dp1[i - 1]
            if((starts1[i - 1] == starts1[i]) ^ (ends1[i - 1] == ends1[i])){
                if(j >= 2) dp1[i][j] = max(dp1[i][j], dp1[i - 1][j - 2] + maxes1[i]);
            }
            else if((starts1[i - 1] == starts1[i]) && (ends1[i - 1] == ends1[i])){
                dp1[i][j] = max(dp1[i][j], dp1[i - 1][j] + maxes1[i]);
            }
            else{
                if(j >= 4) dp1[i][j] = max(dp1[i][j], dp1[i - 1][j - 4] + maxes1[i]);
            }

            // From dp2[i - 1]
            if((starts1[i - 1] == starts1[i]) ^ (ends1[i - 1] == ends1[i])){
                if(j >= 2) dp1[i][j] = max(dp1[i][j], dp2[i - 1][j - 2] + maxes1[i]);
            }
            else if((starts1[i - 1] == starts1[i]) && (ends1[i - 1] == ends1[i])){
                dp1[i][j] = max(dp1[i][j], dp2[i - 1][j] + maxes1[i]);
            }
            else{
                if(j >= 4) dp1[i][j] = max(dp1[i][j], dp2[i - 1][j - 4] + maxes1[i]);
            }



            // dp2
            // From dp0[i - 1]
            if (j >= 10) {
                dp2[i][j] = max(dp2[i][j], dp0[i-1][j-10] + maxes1[i] + maxes2[i]);
            }

            // From dp1[i - 1]
            if(starts1[i - 1] == starts1[i] && ends1[i - 1] == ends1[i]){
                if(j >= 6) dp2[i][j] = max(dp2[i][j], dp1[i - 1][j - 6] + maxes1[i] + maxes2[i]);
            }
            else{
                if(j >= 8) dp2[i][j] = max(dp2[i][j], dp1[i - 1][j - 8] + maxes1[i] + maxes2[i]);
            }

            // From dp2[i - 1]
            if(starts1[i - 1] == starts1[i] && ends1[i - 1] == ends1[i]){
                if(j >= 6) dp2[i][j] = max(dp2[i][j], dp2[i - 1][j - 6] + maxes1[i] + maxes2[i]);
            }
            else{
                if(j >= 8) dp2[i][j] = max(dp2[i][j], dp2[i - 1][j - 8] + maxes1[i] + maxes2[i]);
            }
        }
    }

    int optimum = 0;
    int used = -1;
    for (int i=1; i<=1000; i++){
        int temp = max({dp0[divisions - 1][i], dp1[divisions - 1][i], dp2[divisions - 1][i]});
        if (temp >= optimum){
            optimum = temp;
            used = i;
        }
    }

    return {optimum,used};
}

vector<int> getSubarrays(int rows, int used) {
    int divisions = 10000 / rows;
    int i = -1, j = divisions - 1, k = used;
    vector<int> selected(divisions);
    int mx = max({dp0[j][k], dp1[j][k], dp2[j][k]});
    if(mx == dp0[j][k]) i = 0;
    else if(mx == dp1[j][k]) i = 1;
    else i = 2;
    
    while (j >= 0) {
        int which = i;
        
        selected[j] = i;
        if(j == 0){
            break;
        }

        if (which == 0) {
            int m = dp0[j][k];
            if (m == dp0[j - 1][k]) {
                i = 0;
            } else if (m == dp1[j - 1][k]) {
                i = 1;
            } else {
                i = 2;
            }
            j -= 1;
            continue;
        }
        
        else if (which == 1){
            int m = dp1[j][k];

            // From dp0[j - 1] ???
            if (k >= 6 && m == dp0[j - 1][k - 6] + maxes1[j]){
                i = 0; j -= 1; k -= 6;
                continue;
            }
            // From dp1[j - 1] ???
            if ((starts1[j - 1] == starts1[j]) ^ (ends1[j - 1] == ends1[j])) {
                if (k >= 2 && m == dp1[j - 1][k - 2] + maxes1[j]){
                    i = 1; j -= 1; k -= 2;
                    continue;
                }
            }
            else if (starts1[j - 1] == starts1[j] && ends1[j - 1] == ends1[j]) {
                if(m == dp1[j - 1][k] + maxes1[j]){
                    i = 1; j -= 1;
                    continue;
                }
            }
            else {
                if (k >= 4 && m == dp1[j - 1][k - 4] + maxes1[j]){
                    i = 1; j -= 1; k -= 4;
                    continue;
                }                
            }
            // From dp2[j - 1] ???
            if ((starts1[j - 1] == starts1[j]) ^ (ends1[j - 1] == ends1[j])) {
                if (k >= 2 && m == dp2[j - 1][k - 2] + maxes1[j]){
                    i = 2; j -= 1; k -= 2;
                    continue;
                }
            }
            else if (starts1[j - 1] == starts1[j] && ends1[j - 1] == ends1[j]) {
                if(m == dp2[j - 1][k] + maxes1[j]){
                    i = 2; j -= 1;
                    continue;
                }
            }
            else {
                if (k >= 4 && m == dp2[j - 1][k - 4] + maxes1[j]){
                    i = 2; j -= 1; k -= 4;
                    continue;
                }                
            }
        }
        else {
            int m = dp2[j][k];

            // From dp0[j - 1] ???
            if (k >= 10 && m == dp0[j - 1][k - 10] + maxes1[j] + maxes2[j]) {
                i = 0; j -= 1; k -= 10;
                continue;
            }
            // From dp1[j - 1] ???
            if (starts1[j - 1] == starts1[j] && ends1[j - 1] == ends1[j]) {
                if (k >= 6 && m == dp1[j - 1][k - 6] + maxes1[j] + maxes2[j]){
                    i = 1; j -= 1; k -= 6;
                    continue;
                }
            }
            else {
                if (k >= 8 && m == dp1[j - 1][k - 8] + maxes1[j] + maxes2[j]){
                    i = 1; j -= 1; k -= 8;
                    continue;
                }
            }
            // From dp2[j - 1] ???
            if (starts1[j - 1] == starts1[j] && ends1[j - 1] == ends1[j]) {
                if (k >= 6 && m == dp2[j - 1][k - 6] + maxes1[j] + maxes2[j]){
                    i = 2; j -= 1; k -= 6;
                    continue;
                }
            }
            else {
                if (k >= 8 && m == dp2[j - 1][k - 8] + maxes1[j] + maxes2[j]){
                    i = 2; j -= 1; k -= 8;
                    continue;
                }
            }
        }
    }
    return selected;

}

bool point_cmp(POINT a, POINT b){
    return (a.first == b.first && a.second == b.second);
}

edgesStruct getEdges(vector<int> selected, int len){
    
    vector<EDGE> edgeList;
    vector<POINT> vertices;
    vector<EDGE> upEdges, downEdges;

    int n = starts1.size();

    int initial = 0;

    while(selected[initial] == 0){
        initial++;
    }

    int x = initial*len;
    POINT preStart, preEnd;

    if(selected[initial]==1){
        POINT A,B;
        A = {x-0.6,starts1[initial]-0.6};
        B = {x-0.6,ends1[initial]+0.6};
        vpb(A);
        vpb(B);
        epb({A,B});
        preStart = A;
        preEnd = B;
    }
    
    else if(selected[initial]==2){

        POINT A = {x+len-1+0.2,starts2[initial]-0.6}; // right side start of small block
        POINT B = {x+len-1+0.2,ends2[initial]+0.6}; // right side end of small block
        vpb(A);
        vpb(B);
        epb({A, B});

        // extreme are extended by 0.6 and others by 0.4
        POINT C,D,E,F;
        if(starts2[initial] > ends1[initial]){
            C = {x-0.6,starts1[initial]-0.6};
            D = {x-0.6,ends2[initial]+0.6};
            E = {x-0.4,starts2[initial]-0.6};
            F = {x-0.4,ends1[initial]+0.6};            
            vpb(C); 
            vpb(D); 
            vpb(E); 
            vpb(F);
            epb({C, D}); 
            epb({E, F}); 
            epb({B, D}); 
            epb({A, E});
            preStart = C;
            preEnd = F;
        }
        else{
            C = {x-0.6,starts2[initial]-0.6};
            D = {x-0.6,ends1[initial]+0.6};
            E = {x-0.4,starts1[initial]-0.6};
            F = {x-0.4,ends2[initial]+0.6};            
            vpb(C); 
            vpb(D); 
            vpb(E); 
            vpb(F);
            epb({C, D}); 
            epb({E, F}); 
            epb({F, B}); 
            epb({C, A});
            preStart = E;
            preEnd = D;
        }
    }

    int previndex = initial;

    bool iszero = 0;

    for(int i = initial+1;i<n;i++){
        int curr = selected[i];
        if(curr == 0){
            iszero = 1;
            continue;
        }

        int x = i*len;
        int prevx = previndex*len;

        if(iszero){
            iszero = 0;
            // previous large box's rightmost vertices and edges
            POINT X = {prevx+len-1+0.6,starts1[previndex]-0.6};
            POINT Y = {prevx+len-1+0.6,ends1[previndex]+0.4}; //0.2 for pipeline
            
            vpb(X);
            vpb(Y); 
            // previous box remaining 2 edges
            epb({X,Y});
            epb({X, preStart});


            if(curr == 2){

                // right side of small block
                POINT A = {x+len-1+0.2,starts2[i]-0.6};
                POINT B = {x+len-1+0.2,ends2[i]+0.6};
                vpb(A);
                vpb(B);
                epb({A,B});

                POINT C,D,E,F,G,H;

                if(ends1[previndex] > max(ends1[i],ends2[i])){
                    C = {x-0.4,ends1[previndex]+0.6};
                    D = {x-0.6,ends1[previndex]+0.4};
                    E = {x-0.4, max(ends1[i], ends2[i]) + 0.6};
                    G = {x-0.4, min(ends1[i], ends2[i]) + 0.6};
                    F = {x-0.4,max(starts1[i], starts2[i])-0.6};
                    H = {x-0.6,min(starts1[i], starts2[i])-0.6};
                    vpb(C); 
                    vpb(D); 
                    vpb(E); 
                    vpb(F); 
                    vpb(G); 
                    vpb(H);
                    epb({F,G}); 
                    epb({C,E}); 
                    epb({D,H}); 
                    epb({D,Y}); 
                    epb({C,preEnd});

                    upEdges.push_back({preEnd,C});
                    downEdges.push_back({Y,D});

                    if(ends2[i]>ends1[i]){
                        epb({E,B}); 
                        epb({F,A});
                        preEnd = G;
                        preStart = H;
                    }
                    else{
                        epb({G,B}); 
                        epb({H,A});
                        preEnd = E;
                        preStart = F;
                    }

                }
                else if(ends1[previndex] < min(starts1[i],starts2[i])){
                    C = {x-0.4,ends1[previndex]+0.4};
                    D = {x-0.6,ends1[previndex]+0.6};
                    E = {x-0.4,min(starts1[i], starts2[i]) -0.6};
                    F = {x-0.4,min(ends1[i], ends2[i]) +0.6};
                    G = {x-0.4,max(starts1[i], starts2[i]) -0.6};
                    H = {x-0.6,max(ends1[i], ends2[i])+0.6};
                    
                    vpb(C); 
                    vpb(D); 
                    vpb(E); 
                    vpb(F); 
                    vpb(G); 
                    vpb(H);
                    epb({D, preEnd});  
                    epb({Y,C}); 
                    epb({E,C}); 
                    epb({D,H}); 
                    epb({G,F});

                    upEdges.push_back({preEnd,D});
                    downEdges.push_back({Y,C});

                    if(ends2[i]>ends1[i]){
                        epb({B,H}); 
                        epb({A,G});
                        preEnd = F;
                        preStart = E;
                    }
                    else{
                        epb({F,B}); 
                        epb({E,A});
                        preEnd = H;
                        preStart = G;
                    }

                }
                else{
                    C = {x-0.6,max(ends1[i], ends2[i])+0.6};
                    F = {x-0.6,min(starts1[i], starts2[i])-0.6};
                    E = {x-0.4,min(ends1[i], ends2[i])+0.6};
                    D = {x-0.4,max(starts1[i],starts2[i])-0.6};
                    H = {x-0.6,ends1[previndex]+0.6};
                    G = {x-0.6,ends1[previndex]+0.4};
                    vpb(C); 
                    vpb(D); 
                    vpb(E); 
                    vpb(F); 
                    vpb(G); 
                    vpb(H);
                    epb({preEnd,H}); 
                    epb({Y,G}); 
                    epb({H,C}); 
                    epb({F,G}); 
                    epb({D,E}); 

                    upEdges.push_back({preEnd,H});
                    downEdges.push_back({Y,G});

                    if(ends2[i]>ends1[i]){
                        epb({C,B}); epb({D,A});
                        preEnd = E;
                        preStart = F;
                    }
                    else{
                        epb({E,B}); epb({F,A});
                        preEnd = C;
                        preStart = D;
                    }

                }

            }

            if(curr == 1){
                POINT C,D,E,F;
                if(ends1[previndex] > ends1[i]){
                    C = {x-0.6,ends1[previndex]+0.4};
                    D = {x-0.4,ends1[previndex]+0.6};
                    E = {x-0.4,ends1[i]+0.6};
                    F = {x-0.6,starts1[i]-0.6};

                    vpb(C); 
                    vpb(D); 
                    vpb(E); 
                    vpb(F);
                    epb({preEnd,D}); 
                    epb({Y,C}); 
                    epb({D,E}); 
                    epb({C,F});
                    preEnd = E;
                    preStart = F;

                    upEdges.push_back({preEnd,D});
                    downEdges.push_back({Y,C});
                }
                else if(ends1[previndex] < starts1[i]){
                    C = {x-0.4,ends1[previndex]+0.4};
                    D = {x-0.6,ends1[previndex]+0.6};
                    E = {x-0.4,starts1[i]-0.6};
                    F = {x-0.6,ends1[i]+0.6};

                    vpb(C); 
                    vpb(D); 
                    vpb(E); 
                    vpb(F);
                    epb({preEnd,D}); 
                    epb({Y,C}); 
                    epb({C,E}); 
                    epb({D,F});
                    preEnd = F;
                    preStart = E;

                    upEdges.push_back({preEnd,D});
                    downEdges.push_back({Y,C});
                }
                else{
                    C = {x-0.6,ends1[i]+0.6};
                    D = {x-0.6,ends1[previndex]+0.6};
                    E = {x-0.6,ends1[previndex]+0.4};
                    F = {x-0.6,starts1[i]-0.6};
                    
                    vpb(C); 
                    vpb(D); 
                    vpb(E); 
                    vpb(F);
                    epb({preEnd,D}); 
                    epb({Y,E}); 
                    epb({D,C}); 
                    epb({E,F});
                    preEnd = C;
                    preStart = F;

                    upEdges.push_back({preEnd,D});
                    downEdges.push_back({Y,E});
                }
            }

        }
        
        else{
            if(curr==1){
                POINT C,D,E,F;

                if(starts1[previndex] == starts1[i] && (ends1[previndex] == ends1[i])){
                    previndex = i;
                    continue;
                }

                if(starts1[previndex] == starts1[i]){
                    if(ends1[i]>ends1[previndex]){
                        C = {x-0.6, ends1[i]+0.6};
                        D = {x-0.6, ends1[previndex]+0.6};
                    }
                    else{
                        C = {x-0.4, ends1[i]+0.6};
                        D = {x-0.4, ends1[previndex]+0.6};
                    }
                    vpb(C); vpb(D);
                    epb({C,D}); 
                    epb({preEnd,D});
                    preEnd = C;

                }
                else if(ends1[previndex] == ends1[i]){
                    if(starts1[i]>starts1[previndex]){
                        C = {x-0.4, starts1[i]-0.6};
                        D = {x-0.4, starts1[previndex]-0.6};
                    }
                    else{
                        C = {x-0.6, starts1[i]-0.6};
                        D = {x-0.6, starts1[previndex]-0.6};
                    }
                    vpb(C); vpb(D);
                    epb({C,D}); 
                    epb({D, preStart});
                    preStart = C;
                }
                else{
                    if(ends1[i] > ends1[previndex]){
                        C = {x-0.6, ends1[i]+0.6};
                        D = {x-0.6, ends1[previndex]+0.6};

                        if(starts1[i] > starts1[previndex]){
                            E = {x-0.4, starts1[i]-0.6};
                            F = {x-0.4, starts1[previndex]-0.6};
                        }
                        else{
                            E = {x-0.6, starts1[i]-0.6};
                            F = {x-0.6, starts1[previndex]-0.6};
                        }

                        vpb(C); vpb(D); vpb(E); vpb(F);
                        epb({C,D}); 
                        epb({E,F}); 
                        epb({preEnd,D}); 
                        epb({preStart,F});
                        preStart = E;
                        preEnd = C;
                    }
                    else{
                        C = {x-0.4, ends1[previndex]+0.6};
                        D = {x-0.4, ends1[i]+0.6};

                        if(starts1[i] > starts1[previndex]){
                            E = {x-0.4, starts1[previndex]-0.6};
                            F = {x-0.4, starts1[i]-0.6};
                        }
                        else{
                            E = {x-0.6, starts1[previndex]-0.6};
                            F = {x-0.6, starts1[i]-0.6};
                        }

                        vpb(C); vpb(D); vpb(E); vpb(F);
                        epb({C,D}); 
                        epb({E,F}); 
                        epb({preEnd,C}); 
                        epb({preStart,E});
                        preStart = F;
                        preEnd = D;
                    }

                }

            }

            else{
                // curr = 2;
                POINT A,B,C,D,E,F,G,H;
                A = {x+len-1+0.2, starts2[i]-0.6};
                B = {x+len-1+0.2, ends2[i]+0.6};

                C = {x-0.4, max(starts1[i], starts2[i]) -0.6};
                D = {x-0.4, min(ends1[i], ends2[i]) +0.6};

                vpb(A); 
                vpb(B); 
                vpb(C); 
                vpb(D);
                epb({A,B}); 
                epb({C,D});

                double newStart = min(starts1[i], starts2[i]);
                double newEnd = max(ends1[i], ends2[i]);

                if(newEnd > ends1[previndex]){
                    E = {x-0.6, newEnd+0.6};
                    F = {x-0.6, ends1[previndex]+0.6};

                    if(newStart > starts1[previndex]){
                        G = {x-0.4, newStart-0.6};
                        H = {x-0.4, starts1[previndex]-0.6};
                    }
                    else{
                        G = {x-0.6, newStart-0.6};
                        H = {x-0.6, starts1[previndex]-0.6};
                    }
                    vpb(E); 
                    vpb(F); 
                    vpb(G); 
                    vpb(H);
                    epb({preEnd,F}); 
                    epb({preStart,H});  
                    epb({G,H}); 
                    epb({E,F});
                    
                    if(ends2[i]>ends1[i]){
                        epb({E,B}); 
                        epb({C,A});
                        preStart = G;
                        preEnd = D;
                    }
                    else{
                        epb({D,B}); 
                        epb({G,A});
                        preStart = C;
                        preEnd = E;
                    }
                }
                else{
                    E = {x-0.4, ends1[previndex]+0.6};
                    F = {x-0.4, newEnd+0.6};

                    if(newStart > starts1[previndex]){
                        G = {x-0.4, newStart-0.6};
                        H = {x-0.4, starts1[previndex]-0.6};
                    }
                    else{
                        G = {x-0.6, newStart-0.6};
                        H = {x-0.6, starts1[previndex]-0.6};
                    }

                    vpb(E); 
                    vpb(F); 
                    vpb(G); 
                    vpb(H);
                    epb({preEnd,E}); 
                    epb({preStart,H});  
                    epb({G,H}); 
                    epb({E,F});

                    if(ends2[i]>ends1[i]){
                        epb({F,B}); 
                        epb({C,A});
                        preStart = G;
                        preEnd = D;
                    }
                    else{
                        epb({D,B}); 
                        epb({G,A});
                        preStart = C;
                        preEnd = F;
                    }


                }

            }
        }

        previndex = i;
    }
    POINT lastS, lastE;
    int lastX = len*previndex + len - 1 + 0.6;
    lastS = {lastX, starts1[previndex]-0.6};
    lastE = {lastX, ends1[previndex]+0.6};
    epb({lastE, lastS});
    epb({preEnd, lastE});
    epb({preStart, lastS});

    vpb(lastE); vpb(lastS);

    return {edgeList , upEdges , downEdges , vertices};
}

void extend_edges(edgesStruct & es ){

    set<EDGE> anslist;
    for (auto &&i : es.edgeList)
    {
        anslist.insert(i);
    }
    
    map<POINT , vector<EDGE>> adj;
    for (auto &&i : es.edgeList)
    {
        adj[i.first].push_back(i);
        adj[i.second].push_back(i);
    }

    for (auto &&i : es.upEdges)
    {
        float cury = i.first.second;
        float x1 = i.first.first;
        float x2 = i.second.first;
        if( x1 == 0 or x2 == 0 or floor(cury) == 0) continue;
        int max_sum  = 0;
        float maxy = -1;
        for (int j = ceil(cury) ; j < 10002 ; j++)
        {
            int this_sum = getSum(x1 , floor(cury) , x2 , j);
            if(this_sum > max_sum){
                max_sum = this_sum;
                maxy = j;
            }
        }

        if(maxy != -1){
            for (auto &&j : adj[{x1, cury}])
            {
                anslist.erase(j);
                POINT p2 = ( point_cmp(j.first , {x1, cury}) ? j.second : j.first);
                if( !point_cmp(p2, {x2 , cury}) ) anslist.insert({{x1, maxy}, p2});
            }
            for (auto &&j : adj[{x2, cury}])
            {
                anslist.erase(j);
                POINT p2 = ( point_cmp(j.first , {x2, cury}) ? j.second : j.first);
                if( !point_cmp(p2, {x1 , cury}) ) anslist.insert({{x2, maxy}, p2});
            }
            anslist.insert({{x1, maxy}, {x2, maxy}});
        }
    }
    
    for (auto &&i : es.downEdges)
    {
        float cury = i.first.second;
        float x1 = i.first.first;
        float x2 = i.second.first;
        if( x1 == 0 or x2 == 0 or floor(cury) == 0) continue;
        int max_sum  = 0;
        float maxy = -1;
        for (int j = floor(cury) ; j > 0 ; j--)
        {
            int this_sum = getSum(x1 , floor(cury) , x2 , j);
            if(this_sum > max_sum){
                max_sum = this_sum;
                maxy = j;
            }
        }

        if(maxy != -1){
            for (auto &&j : adj[{x1, cury}])
            {
                anslist.erase(j);
                POINT p2 = ( point_cmp(j.first , {x1, cury}) ? j.second : j.first);
                if( !point_cmp(p2, {x2 , cury}) ) anslist.insert({{x1, maxy}, p2});
            }
            for (auto &&j : adj[{x2, cury}])
            {
                anslist.erase(j);
                POINT p2 = ( point_cmp(j.first , {x2, cury}) ? j.second : j.first);
                if( !point_cmp(p2, {x1 , cury}) ) anslist.insert({{x2, maxy}, p2});
            }
            anslist.insert({{x1, maxy}, {x2, maxy}});
        }
    }

    vector<EDGE> final_ans;
    for (auto &&i : anslist){
        final_ans.push_back(i);
    }
    es.edgeList = final_ans;
}


// END OF CODE

