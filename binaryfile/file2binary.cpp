#include <iostream>  
#include <fstream>  
#include <vector>  
#include <string>  
#include <iomanip>  

const size_t BUFFER_SIZE = 1024000;
using namespace std;

#define DATA_TYPE double  
void show(vector<double>& x, const char* msg)
{
    std::cout << std::setprecision(15);
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
        std::cerr << "Usage: " << argv[0] << " <input_file_path> <output_file_path>" << std::endl;
        return 1;
    }

    std::string finname = argv[1];
    std::string foutname = argv[2];

    std::ifstream fin(finname, std::ios::binary);
    std::ofstream fout(foutname, std::ios::binary | std::ios::out);

    if (!fin || !fout) {
        std::cerr << "Failed to open file(s)." << std::endl;
        return 1;
    }

    std::vector<DATA_TYPE> buffer(BUFFER_SIZE);
    DATA_TYPE value;
    size_t count = 0;
    size_t num_ts = 0;

    while (fin >> value) {
        num_ts++;
        buffer[count++] = value;

        if (count == BUFFER_SIZE) {
            cout << ".";
            fout.write(reinterpret_cast<char*>(buffer.data()), count * sizeof(DATA_TYPE));
            count = 0;
        }
    }

    if (count > 0) {
        fout.write(reinterpret_cast<char*>(buffer.data()), count * sizeof(DATA_TYPE));
    }

    fin.close();
    fout.close();
    cout << "Number of time series values read: " << num_ts << endl;

    return 0;
}