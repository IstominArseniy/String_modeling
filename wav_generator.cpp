#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

const int sampleRate = 44100;
const int bitDepth = 16;

class SineOscillator {
    float frequency, amplitude, angle = 0.0f, offset = 0.0f;
public:
    SineOscillator(float freq, float amp) : frequency(freq), amplitude(amp) {
        offset = 2 * M_PI * frequency / sampleRate;
    }
    float process() {
        float sample = amplitude * sin(angle);
        angle += offset;
        return sample;
        // Asin(2pif/sr)
    }
};

vector<float> generateRawSineWave(float frequency, float amplitude, float duration){
/*
 frequency in hertz, amplitude should be from 0 to 1, duration in sec
 */
    SineOscillator sineOscillator(frequency,amplitude);
    vector<float> rawData;
    for(int i = 0; i < sampleRate * duration; ++i) {
        float sample = sineOscillator.process();
        rawData.push_back(sample);
    }
    return rawData;
}

void writeToFile(ofstream &file, int value, int size) {
    file.write(reinterpret_cast<const char*> (&value), size);
}

void writeWAV(vector<float> rawData, string filename){
    /*
     rawData should contain floats from -1 to 1
     this function work good for every sampleRate and every bitDepth(%8 = 0) - I don't test it and don't promise it, but it should)
     */
    ofstream audioFile;
    audioFile.open(filename + ".wav", ios::binary);
    //Header chunk
    audioFile << "RIFF";
    audioFile << "----";
    audioFile << "WAVE";
    // Format chunk
    audioFile << "fmt ";
    writeToFile(audioFile, 16, 4); // Size
    writeToFile(audioFile, 1, 2); // Compression code
    writeToFile(audioFile, 1, 2); // Number of channels
    writeToFile(audioFile, sampleRate, 4); // Sample rate
    writeToFile(audioFile, sampleRate * bitDepth / 8, 4 ); // Byte rate
    writeToFile(audioFile, bitDepth / 8, 2); // Block align
    writeToFile(audioFile, bitDepth, 2); // Bit depth
    //Data chunk
    audioFile << "data";
    audioFile << "----";
    //Writing rawData
    int preAudioPosition = audioFile.tellp();

    auto maxAmplitude = pow(2, bitDepth - 1) - 1;
    int dataSize = rawData.size();
    for(int i = 0; i < dataSize; ++i) {
        int intSample = static_cast<int> (rawData[i] * maxAmplitude);
        writeToFile(audioFile, intSample, bitDepth / 8);
    }
    int postAudioPosition = audioFile.tellp();

    audioFile.seekp(preAudioPosition - 4);
    writeToFile(audioFile, postAudioPosition - preAudioPosition, 4);

    audioFile.seekp(4, ios::beg);
    writeToFile(audioFile, postAudioPosition - 8, 4);

    audioFile.close();
}
int main() {
    vector<float> data;
    data = generateRawSineWave(440, 0.5, 2);
    writeWAV(data, "mySine");
    return 0;
}