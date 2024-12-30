#include <iostream>  
#include <fstream>  
#include <vector>  
#include <string>  

const size_t DATA_SIZE = sizeof(double);

using namespace std;

bool readBinaryData(const string& filename, vector<double>& data) {
    ifstream finn(filename, ios::binary | ios::in);
    if (!finn) {
        cerr << "Failed to open file for reading: " << filename << endl;
        return false;
    }

    finn.seekg(0, ios::end);
    std::streampos fileSize = finn.tellg();
    finn.seekg(0, ios::beg);

    size_t numElements = fileSize / DATA_SIZE;
    data.resize(numElements);

    if (!finn.read(reinterpret_cast<char*>(data.data()), fileSize)) {
        cerr << "Failed to read data from file: " << filename << endl;
        finn.close();
        return false;
    }

    finn.close();
    return true;
}

bool writeBinaryData(const string& filename, const vector<double>& data) {
    ofstream fout(filename);
    if (!fout) {
        cerr << "Failed to open file for writing: " << filename << endl;
        return false;
    }

    for (auto num : data) 
    {
        fout << num << " ";
    }

    fout.close();
    return true;
}

void show(const vector<double>& x, const char* msg) {
    cout << msg << endl;
    cout << "Number of elements: " << x.size() << endl;
    for (size_t i = 0; i < x.size() && i < 10; i++)
        cout << x.at(i) << ", ";
    if (x.size() > 10)
        cout << "...";
    cout << endl;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input_file_path> <output_file_path>" << endl;
        return 1;
    }

    string inputFilename = argv[1];
    string outputFilename = argv[2];
    vector<double> data;

    if (readBinaryData(inputFilename, data)) {
        show(data, "Vector read from binary file.");
        if (writeBinaryData(outputFilename, data)) {
            cout << "Data successfully written to output file: " << outputFilename << endl;
        }
        else {
            cerr << "Error writing data to output file." << endl;
        }
    }
    else {
        cerr << "Error reading data from input file." << endl;
    }

    return 0;
}