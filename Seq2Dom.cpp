#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <random>
#include <cstdlib>
#include <fstream>
#include <regex>
#include <unordered_set>
#include <future>
#include <sys/stat.h>

int maxKeyDomLength = 6; //independent variable.


using namespace std;

map<char, char> complement;



string stringToBinary(const string &str) {
    string binaryString = "";
    for (char c : str) {
        // Convert each character to a binary representation of ASCII
        binaryString += std::bitset<8>(c).to_string();
    }
    return binaryString;
}


//assume even-number length.
string binary2Base(string binarySequence){
    string baseSequence = "";
    for(int i = 0; i < binarySequence.length(); i+=2){
        string bit = binarySequence.substr(i, 2);
        if(bit == "00") baseSequence += "A";
        else if(bit == "01") baseSequence += "C";
        else if(bit == "10") baseSequence += "G";
        else baseSequence += "T";
    }
    return baseSequence;
}



int countDistinctReactions() {
    string filePath = "peppercornSim/enum.pil";
    ifstream file(filePath);
    string line;
    // This regex matches both "StrandX + StrandXKEY -> ..." and "StrandXKEY + StrandX -> ..."
    regex reactionPattern(R"(\b(Strand\w+) \+ (Strand\w+KEY) -> \w+|\b(Strand\w+KEY) \+ (Strand\w+) -> \w+)");
    unordered_set<string> uniqueReactions;

    if (!file.is_open()) {
        cerr << "Failed to open file." << endl;
        return -1;
    }

    // Skip lines until "# Condensed reactions"
    bool startCount = false;
    while (getline(file, line)) {
        if (line.find("# Condensed reactions") != string::npos) {
            startCount = true;
            continue;
        }
        if (startCount) {
            smatch matches;
            if (regex_search(line, matches, reactionPattern)) {
                if (matches[1].matched && matches[2].matched) {
                    uniqueReactions.insert(matches[1].str() + " + " + matches[2].str());
                } else if (matches[3].matched && matches[4].matched) {
                    uniqueReactions.insert(matches[3].str() + " + " + matches[4].str());
                }
            }
        }
    }

    file.close();
    return uniqueReactions.size();
}


string generateRandomDomainName(size_t minLength, size_t maxLength) {
        // likely repeat (according to birthday paradox estimate) of sqrt(2*characterSetLength^numberOfCharactersPerString*ln(2))

    const char charset[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    const size_t maxIndex = sizeof(charset) - 2; // -2 to account for the null terminator

    random_device rd;
    mt19937 generator(rd());
    uniform_int_distribution<size_t> lengthDistribution(minLength, maxLength);
    uniform_int_distribution<size_t> charDistribution(0, maxIndex);

    size_t randomLength = lengthDistribution(generator);
    string randomString;
    randomString.reserve(randomLength); // Reserve to optimize memory allocation

    for (size_t i = 0; i < randomLength; ++i) {
        randomString += charset[charDistribution(generator)];
    }

    return randomString;
}




struct Overlap{
    int startI = -1, finishI = -1, startJ = -1, finishJ = -1, strandA = -1, strandB = -1;

    bool operator<(const Overlap& other) const {
        int length = finishI-startI, otherLength = other.finishI-other.startI;
        if(length == otherLength){
            return strandA < other.strandA;
        }
        return length > otherLength;
        
    }
    int length(){
        return finishI-startI;
    }

};
// handle loose characters later.

struct Domain{
    string name;
    int start, end;
    bool complement;
    int strand;
    bool operator<(const Domain& other) const {
        
        if(strand == other.strand){
            return start < other.start;
        }
        return strand < other.strand;
        
    }
    int length(){ 
        return end-start;
    }

};

struct Strand{
    string name;
    string structure;
    int length = 0;
    bool isKey;

    // note: doesn't reset isKey (intentional for the purposes of the key construction).
    void clear(){
        name = "";
        structure = "";
        length = 0;
    }

};


pair<pair<int, int>, pair<int, int> > longestCommonComplementIndices(const string& s1, const string& s2) { 
    int m = s1.length(), n = s2.length();
    vector<vector<int> > dp(m + 1, vector<int>(n + 1, 0));
    int maxLength = 0;
    int endIdx1 = -1;  // End index in s1
    int endIdx2 = -1;  // End index in s2

    // Fill the DP table
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (s1[i - 1] == complement[s2[j - 1]]) {
                dp[i][j] = dp[i - 1][j - 1] + 1;
                if (dp[i][j] > maxLength) {
                    maxLength = dp[i][j];
                    endIdx1 = i - 1;  // Update the end index for s1
                    endIdx2 = j - 1;  // Update the end index for s2
                }
            } else {
                dp[i][j] = 0;
            }
        }
    }

    if (maxLength == 0) {
        return {{-1, -1}, {-1, -1}};  // No common substring found
    }
    // Calculate the start indices and return indices for both strings
    int startIdx1 = endIdx1 - maxLength + 1;
    int startIdx2 = endIdx2 - maxLength + 1;
    return {{startIdx1, endIdx1}, {startIdx2, endIdx2}};
}


int reportPerformance() {
    string filePath = "peppercornSim/enum.pil";
    ifstream file(filePath);
    string line;
    // Adjusted regex to match both StrandX + StrandXKey and StrandXKey + StrandX
    regex reactionPattern(R"(Strand(\w+) \+ Strand\1KEY -> \w+|Strand(\w+)KEY \+ Strand\2 -> \w+)");
    unordered_set<std::string> uniqueReactions;

    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return -1;
    }

    // Skip lines until "# Condensed reactions"
    bool startCount = false;
    while (getline(file, line)) {
        if (line.find("# Condensed reactions") != std::string::npos) {
            startCount = true;
            continue;
        }
        if (startCount) {
            smatch matches;
            while (regex_search(line, matches, reactionPattern)) {
                if (matches[1].matched) {  // Checks the first group
                    uniqueReactions.insert("Strand" + matches[1].str() + " + Strand" + matches[1].str() + "KEY");
                } else if (matches[2].matched) {  // Checks the second group
                    uniqueReactions.insert("Strand" + matches[2].str() + "KEY + Strand" + matches[2].str());
                }
                line = matches.suffix().str();  // Move to the next potential match
            }
        }
    }

    file.close();
    return uniqueReactions.size();
}


double calculateAverageRate() {
    string filePath = "peppercornSim/enum.pil";
    ifstream file(filePath);
    string line;
    regex ratePattern(R"(\[condensed\s+=\s+([0-9.]+) /nM/s \])");  // Regex to capture the rate value
    double sum = 0.0;
    int count = 0;

    if (!file.is_open()) {
        cerr << "Failed to open file." << endl;
        return -1;  // Return -1 to indicate error
    }

    bool startCount = false;
    while (getline(file, line)) {
        if (line.find("# Condensed reactions") != string::npos) {
            startCount = true;
            continue;
        }
        if (startCount) {
            smatch matches;
            if (regex_search(line, matches, ratePattern) && matches.size() > 1) {
                sum += stod(matches[1]);  // Convert string to double and add to sum
                count++;
            }
        }
    }

    file.close();
    
    if (count == 0) return 0;  // Prevent division by zero
    
    return sum / count;  // Calculate average
}


void runPeppercorn(){
    system("peppercorn -o peppercornSim/enum.pil peppercornSim/main.pil -c  -L 7 --max-complex-size 2  --ignore-branch-4way ");
}



bool fileExists(const std::string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

void writeToCSV(const std::string& filePath, int r, double successfulReactionPercent, 
                int numberOfFiles, double avgRate, double totalAverageTime, double keyToFileRatio, 
                double avgFileLength, int maxKeyDomLength) {
    std::ofstream outputFile;
    bool exists = fileExists(filePath);

    // Open file for appending, create if it does not exist
    outputFile.open(filePath, std::ios::out | std::ios::app);

    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        return;
    }

    // If the file is new or empty, write the headers
    if (!exists) {
        outputFile << "Maximum Key Length,Successful Reaction %,Number Of Files,"
                      "Average Time of Reaction,Average Time of Full Consumption,Average Key-to-file Ratio,"
                      "Average File Length,Successful Reactions\n";
    }

    // Write data
    outputFile << maxKeyDomLength << ","
                << successfulReactionPercent << ","
                << numberOfFiles << ","
                << avgRate << ","
                << totalAverageTime << ","
                << keyToFileRatio << ","
                << avgFileLength << ","
                << r << "\n";

    outputFile.close();
}




int main(){
    freopen("textFiles/input.txt", "r", stdin);
    freopen("peppercornSim/main.pil", "w", stdout);
    vector<Domain> domains;

    complement['A'] = 'T';
    complement['T'] = 'A';
    complement['C'] = 'G';
    complement['G'] = 'C';

    int numOfFiles; cin >> numOfFiles;
    vector<string> binaryFiles (numOfFiles);
    //problematic
    string stringFile;
    for(int i = 0; i < numOfFiles; i++){
        cin >> stringFile;
        binaryFiles[i] = stringToBinary(stringFile);
    }
    for(int i = 0; i < binaryFiles.size(); i++){
        binaryFiles[i] = binary2Base(binaryFiles[i]);
    }

    for(int i = 0; i < binaryFiles.size(); i++){
        vector<Overlap> overlaps (binaryFiles.size());
        for(int j = 0; j < binaryFiles.size(); j++){

            string fileA = binaryFiles[i];
            string fileB = binaryFiles[j];
            pair<pair<int, int>, pair<int, int> > lcsInfo = longestCommonComplementIndices(fileA, fileB);
            int lcsLength = lcsInfo.first.second - lcsInfo.first.first;
            int currentMaxLength = overlaps[i].length();
            if(lcsLength > currentMaxLength){
                overlaps[i].startI = lcsInfo.first.first;
                overlaps[i].finishI = lcsInfo.first.second;
                overlaps[i].startJ = lcsInfo.second.first;
                overlaps[i].finishJ = lcsInfo.second.second;
                overlaps[i].strandA = i;
                overlaps[i].strandB = j;

            }
        }
        sort(overlaps.begin(), overlaps.end());
        Overlap strongestAvailable = overlaps[0];
        string randomDomainName = generateRandomDomainName(3,5);
        if(strongestAvailable.startI != -1){
            Domain dom; 
            dom.start = strongestAvailable.startI;
            dom.end = strongestAvailable.finishI;
            dom.name = randomDomainName;
            dom.strand = strongestAvailable.strandA;
            dom.complement = false;
            domains.push_back(dom);
        }
        if(strongestAvailable.startJ != -1){
            Domain domAsterisk;
            domAsterisk.start = strongestAvailable.startJ; //patched.
            domAsterisk.end = strongestAvailable.finishJ; // patched.
            domAsterisk.name = randomDomainName + "*";
            domAsterisk.strand = strongestAvailable.strandB;
            domAsterisk.complement = true;
            domains.push_back(domAsterisk);
        }

        for(int i = strongestAvailable.startI; i <= strongestAvailable.finishI; i++){ // patched.
            if(i!=-1) binaryFiles[strongestAvailable.strandA][i] = 'X';
        }
        for(int i = strongestAvailable.startJ; i <= strongestAvailable.finishJ; i++){ //patched.
            if(i!=-1) binaryFiles[strongestAvailable.strandB][i] = 'X';
        }
    }

    // remaining nucleotides.
    for(int k = 0; k < binaryFiles.size(); k++){
        string &file = binaryFiles[k];
        int start = 0, end = -1;
        for(int i = 0; i < file.length(); i++){
            
            if(file[i]=='X' || i == file.length()-1){
                if(file[i] != 'X' && i == file.length()-1) end++;
                if(end >= start){
                    
                    Domain looseDomain;
                    looseDomain.complement = false;
                    looseDomain.start = start;
                    looseDomain.end = end;
                    looseDomain.strand = k;
                    looseDomain.name = generateRandomDomainName(3,5);

                    domains.push_back(looseDomain);
                }
                start = i+1;
                end = i;
            }
            if(file[i] != 'X'){
                end++;
                file[i] = 'X';
            }
        }
    }

    sort(domains.begin(), domains.end());
   
    
    vector<Strand> strandS;
    int currentStrandIndex = 0;
    Strand currentStranD, keyStranD;
    currentStranD.isKey = false; keyStranD.isKey = true;
    double avgFileLength = 0;
    double keyToFileRatio = 0;
    for(int i = 0; i < domains.size(); i++){
        auto dom = domains[i];
        if(!(dom.complement)) cout << "length "  << dom.name << " = " << dom.length() << endl;
        if(dom.strand > currentStrandIndex){
            // could define naming here.
            strandS.push_back(currentStranD);
            strandS.push_back(keyStranD);
            avgFileLength += currentStranD.length;
            keyToFileRatio += (double(keyStranD.length) / double(currentStranD.length));

            keyStranD.clear();
            currentStranD.clear();
            currentStrandIndex = dom.strand;
        }
        currentStranD.structure += dom.name + " ";
        currentStranD.length++;
        if(currentStranD.length < maxKeyDomLength){
            
            if(dom.complement){
                keyStranD.structure += dom.name.substr(0, dom.name.length()-1);
                keyStranD.structure += " ";
                keyStranD.length++;
            }
            else{
                keyStranD.structure += dom.name;
                keyStranD.structure += "* ";
                keyStranD.length++;
            }
        }
        if(i == domains.size()-1){
            keyToFileRatio += (double(keyStranD.length) / double(currentStranD.length));
            avgFileLength += currentStranD.length;
            strandS.push_back(currentStranD);
            strandS.push_back(keyStranD);
        }


    }

    keyToFileRatio /= (double(strandS.size())/2);
    avgFileLength /= (strandS.size()/2);
    

    for(int i = 0; i < strandS.size(); i++){
        auto strand = strandS[i];
        string s_name = "Strand";
        
        if(strand.isKey){
            s_name.push_back('A'+i-1);
            s_name += "KEY";
        }
        else{
            s_name.push_back('A'+i);
        }
        
        cout << s_name << " = " << strand.structure << "@initial 100nM" << endl;
    }
    // Pending: define keys.


    

    std::future<void> answer = std::async(std::launch::async, runPeppercorn);

    int r = reportPerformance();
    double avgRate = calculateAverageRate();
    double successfulReactionPercent = double(r) / double(strandS.size()) * 200;
    int numberOfFiles = strandS.size() / 2;
    double totalAverageTime = 100 / avgRate;

    cout << "----------------------";
    cout << "\n" << "Succesful Reactions: " << r << "\n";
    cout << "\n" << "Succesful Reaction %: " << double(r) / double(strandS.size()) * 2 << "\n";
    cout << "\n" << "Number Of Files: " << strandS.size() / 2 << "\n";
    cout << "\n" << "Average Time of Reaction: " << avgRate << "\n";
    cout << "\n" << "Average Key-to-file Ratio: " << keyToFileRatio << "\n";
    cout << "\n" << "Average File Length: " << avgFileLength << "\n";
     cout << "\n" << "Average Total Consumption Time: " << totalAverageTime << "\n";
    cout << "\nMaximum Key Length: " << maxKeyDomLength << "\n";

    string filePath = "output.csv";
    writeToCSV(filePath, r, successfulReactionPercent, numberOfFiles, avgRate, totalAverageTime,
               keyToFileRatio, avgFileLength, maxKeyDomLength);


}
