/*
 * Copyright 2020, Yun (Leo) Zhang <imzhangyun@gmail.com>
 *
 * This file is part of HISAT-3N.
 *
 * HISAT-3N is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT-3N is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT-3N.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <chrono>
#include <cstddef>
#include <iostream>
#include <getopt.h>
#include <stdexcept>
#include <thread>
#include "position_3n_table.h"
#include "utility_3n_table.h"

using namespace std;

string alignmentFileName;
bool standardInMode = false;
string refFileName;
string outputFileName;
bool uniqueOnly = false;
bool multipleOnly = false;
bool CG_only = false;
int nThreads = 1;
long long int loadingBlockSize = 1000000;
char convertFrom = '0';
char convertTo = '0';
char convertFromComplement;
char convertToComplement;
bool addedChrName = false;
bool removedChrName = false;


Positions* positions;

bool fileExist (string& filename) {
    ifstream file(filename);
    return file.good();
}

enum {
    ARG_ADDED_CHRNAME = 256,
    ARG_REMOVED_CHRNAME
};

static const char *short_options = "s:r:t:b:umcp:h";
static struct option long_options[] {
                {"alignments",  required_argument, 0, 'a'},
                {"ref",  required_argument, 0, 'r'},
                {"output-name", required_argument, 0, 'o'},
                {"base-change", required_argument, 0, 'b'},
                {"unique-only", no_argument, 0, 'u'},
                {"multiple-only", no_argument, 0, 'm'},
                {"CG-only", no_argument, 0, 'c'},
                {"threads", required_argument, 0, 'p'},
                {"added-chrname", no_argument, 0, ARG_ADDED_CHRNAME },
                {"removed-chrname", no_argument, 0, ARG_REMOVED_CHRNAME },
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
        };

static void printHelp(ostream& out) {
    out << "hisat-3n-table developed by Yun (Leo) Zhang" << endl;
    out << "Usage:" << endl
        << "hisat-3n-table [options]* --alignments <alignmentFile> --ref <refFile> --output-name <outputFile> --base-change <char1,char2>" << endl
        << "  <alignmentFile>           SORTED SAM filename. Please enter '-' for standard input." << endl
        << "  <refFile>                 reference file (should be FASTA format)." << endl
        << "  <outputFile>              file name to save the 3n table (tsv format). By default, alignments are written to the “standard out” or “stdout” filehandle (i.e. the console)." << endl
        << "  <chr1,chr2>               the char1 is the nucleotide converted from, the char2 is the nucleotide converted to." << endl;
    out << "Options (defaults in parentheses):" << endl
        << " Input:" << endl
        << "  -u/--unique-only          only count the base which is in unique mapped reads." << endl
        << "  -m/--multiple-only        only count the base which is in multiple mapped reads." << endl
        << "  -c/--CG-only              only count CG and ignore CH in reference." << endl
        << "  --added-chrname           please add this option if you use --add-chrname during HISAT-3N alignment." << endl
        << "  --removed-chrname         please add this option if you use --remove-chrname during HISAT-3N alignment." << endl
        << "  -p/--threads <int>        number of threads to launch (1)." << endl
        << "  -h/--help                 print this usage message." << endl;
}

static void parseOption(int next_option, const char *optarg) {
    switch (next_option) {
        case 'a': {
            alignmentFileName = optarg;
            if (alignmentFileName == "-") {
                standardInMode = true;
                break;
            }
            if (!fileExist(alignmentFileName)) {
                cerr << "The alignment file is not exist." << endl;
                throw (1);
            }
            break;
        }
        case 'r': {
            refFileName = optarg;
            if (!fileExist(refFileName)) {
                cerr << "reference (FASTA) file is not exist." << endl;
                throw (1);
            }
            break;
        }
        case 'o':
            outputFileName = optarg;
            break;
        case 'b': {
            string arg = optarg;
            if (arg.size() != 3 || arg[1] != ',') {
                cerr << "Error: expected 2 comma-separated "
                     << "arguments to --base-change option (e.g. C,T), got " << arg << endl;
                throw 1;
            }
            convertFrom = toupper(arg.front());
            convertTo = toupper(arg.back());
            break;
        }
        case 'u':{
            uniqueOnly = true;
            break;
        }
        case 'm': {
            multipleOnly = true;
            break;
        }
        case 'c': {
            CG_only = true;
            break;
        }
        case 'h': {
            printHelp(cerr);
            throw 0;
        }
        case 'p': {
            nThreads = stoi(optarg);
            if (nThreads < 1) {
                nThreads = 1;
            }
            break;
        }
        case ARG_ADDED_CHRNAME: {
            addedChrName = true;
            break;
        }
        case ARG_REMOVED_CHRNAME: {
            removedChrName = true;
            break;
        }
        default:
            printHelp(cerr);
            throw 1;
    }
}

static void parseOptions(int argc, const char **argv) {
    int option_index = 0;
    int next_option;
    while (true) {
        next_option = getopt_long(argc, const_cast<char **>(argv), short_options,
                                  long_options, &option_index);
        if (next_option == -1)
            break;
        parseOption(next_option, optarg);
    }

    // check filenames
    if (refFileName.empty() || alignmentFileName.empty()) {
        cerr << "No reference or SAM file specified!" << endl;
        printHelp(cerr);
        throw 1;
    }

    // give a warning for CG-only
    if (CG_only) {
        if (convertFrom != 'C' || convertTo != 'T') {
            cerr << "Warning! You are using CG-only mode. The the --base-change option is set to: C,T" << endl;
            convertFrom = 'C';
            convertTo = 'T';
        }
    }

    // check if --base-change is empty
    if (convertFrom == '0' || convertTo == '0') {
        cerr << "the --base-change argument is required." << endl;
        throw 1;
    }

    if(removedChrName && addedChrName) {
        cerr << "Error: --removed-chrname and --added-chrname cannot be used at the same time" << endl;
        throw 1;
    }

    // set complements
    convertFromComplement = asc2dnacomp[convertFrom];
    convertToComplement = asc2dnacomp[convertTo];
}

/**
 * give a SAM line, extract the chromosome and position information.
 * return true if the SAM line is mapped. return false if SAM line is not maped.
 */
bool getSAMChromosomePos(string* line, string& chr, long long int& pos) {
    int startPosition = 0;
    int endPosition = 0;
    int count = 0;

    while ((endPosition = line->find("\t", startPosition)) != string::npos) {
        if (count == 2) {
            chr = line->substr(startPosition, endPosition - startPosition);
        } else if (count == 3) {
            pos = stoll(line->substr(startPosition, endPosition - startPosition));
            if (chr == "*") {
                return false;
            } else {
                return true;
            }
        }
        startPosition = endPosition + 1;
        count++;
    }
    return false;
}

/*void opeInFile(ifstream& f) {
    if (alignmentFileName == "-") {
        f = cin;
    } else {
        ifstream alignmentFile;
        alignmentFile.open(alignmentFileName, ios_base::in);
        return alignmentFile;
    }
}*/


int hisat_3n_table()
{
    positions = new Positions(refFileName, nThreads, addedChrName, removedChrName);

    // open #nThreads workers
    vector<thread*> workers;
    for (int i = 0; i < nThreads; i++) {
        workers.push_back(new thread(&Positions::append, positions, i));
    }

    // open a output thread
    thread outputThread;
    outputThread = thread(&Positions::outputFunction, positions, outputFileName);

    // main function, initially 2 load loadingBlockSize (2,000,000) bp of reference, set reloadPos to 1 loadingBlockSize, then load SAM data.
    // when the samPos larger than the reloadPos load 1 loadingBlockSize bp of reference.
    // when the samChromosome is different to current chromosome, finish all sam position and output all.
    ifstream inputFile;
    istream *alignmentFile = &cin;
    if (!standardInMode) {
        inputFile.open(alignmentFileName, ios_base::in);
        alignmentFile = &inputFile;
    }

    string* line; // temporary string to get SAM line.
    string samChromosome; // the chromosome name of current SAM line.
    long long int samPos; // the position of current SAM line.
    long long int reloadPos; // the position in reference that we need to reload.
    long long int lastPos = 0; // the position on last SAM line. compare lastPos with samPos to make sure the SAM is sorted.

    while (alignmentFile->good()) {
        positions->getFreeStringPointer(line);
        if (!getline(*alignmentFile, *line)) {
            positions->returnLine(line);
            break;
        }

        if (line->empty() || line->front() == '@') {
            positions->returnLine(line);
            continue;
        }
        // limit the linePool size to save memory
        while(positions->linePool.size() > 1000 * nThreads) {
            this_thread::sleep_for (std::chrono::microseconds(1));
        }
        // if the SAM line is empty or unmapped, get the next SAM line.
        if (!getSAMChromosomePos(line, samChromosome, samPos)) {
            positions->returnLine(line);
            continue;
        }
        // if the samChromosome is different than current positions' chromosome, finish all SAM line.
        // then load a new reference chromosome.
        if (samChromosome != positions->chromosome) {
            // wait all line is processed
            while (!positions->linePool.empty() || positions->outputPositionPool.size() > 100000) {
                this_thread::sleep_for (std::chrono::microseconds(1));
            }
            positions->appendingFinished();
            positions->moveAllToOutput();
            positions->loadNewChromosome(samChromosome);
            reloadPos = loadingBlockSize;
            lastPos = 0;
        }
        // if the samPos is larger than reloadPos, load 1 loadingBlockSize bp in from reference.
        while (samPos > reloadPos) {
            while (!positions->linePool.empty() || positions->outputPositionPool.size() > 100000) {
                this_thread::sleep_for (std::chrono::microseconds(1));
            }
            positions->appendingFinished();
            positions->moveBlockToOutput();
            positions->loadMore();
            reloadPos += loadingBlockSize;
        }
        if (lastPos > samPos) {
            cerr << "The input alignment file is not sorted. Please use sorted SAM file as alignment file." << endl;
            throw 1;
        }
        positions->linePool.push(line);
        lastPos = samPos;
    }
    //}
    if (!standardInMode) {
        inputFile.close();
    }


    // prepare to close everything.

    // make sure linePool is empty
    while (!positions->linePool.empty()) {
        this_thread::sleep_for (std::chrono::microseconds(100));
    }
    // make sure all workers finished their appending work.
    positions->appendingFinished();
    // move all position to outputPool
    positions->moveAllToOutput();
    // wait until outputPool is empty
    while (!positions->outputPositionPool.empty()) {
        this_thread::sleep_for (std::chrono::microseconds(100));
    }
    // stop all thread and clean
    while(positions->freeLinePool.popFront(line)) {
        delete line;
    }
    positions->working = false;
    for (int i = 0; i < nThreads; i++){
        workers[i]->join();
        delete workers[i];
    }
    outputThread.join();
    delete positions;
    return 0;
}

#include "third_party/mio/mio.hpp"
#include "third_party/BS_thread_pool.hpp"
#include "third_party/concurrentqueue.h"

#include <span>
#include <ranges>
#include <string_view>

class LinesIter: public std::iterator<std::input_iterator_tag, string_view, ptrdiff_t, string_view, string_view> {
public:
	explicit LinesIter(string_view content): content(content) {
        parse_line();        
	}

	LinesIter& operator++() {
        if (line.end() == content.end()) {
            content = content.substr(line.size());
        } else {
            content = content.substr(line.size() + 1);
            parse_line();
        }
        return *this;
	}

    LinesIter operator++(int) {
        auto ret = *this;
        ++(*this);
        return ret;
    }

    reference operator*() const {
        return line;
    }

    bool operator==(const LinesIter &o) const {
        return content.begin() == o.content.begin() && content.end() == o.content.end();
    }

    bool operator!=(const LinesIter &o) const {
        return !(*this == o);
    }

    bool empty() const {
        return content.empty();
    }
    size_t size() const {
        return content.size();
    }

private:
	void parse_line() {
        for (size_t i = 0; i < content.size(); i++) {
            if (content[i] == '\n') {
                line = string_view(content.begin(), content.begin() + i);
                return;
            }
        }
        line = content;
	}

	string_view content;
	string_view line;
};

class ChromosomeDB {
public:
	using Range = pair<size_t, size_t>;
	ChromosomeDB(string inputRefFileName, bool inputAddedChrName, bool inputRemovedChrName): addedChrName(inputAddedChrName), removedChrName(inputRemovedChrName) {
		std::error_code error;
		refFileMap.map(inputRefFileName, error);
		refFile_ = std::span(refFileMap.data(), refFileMap.size());
		assert(!error);
		LoadChromosomeNamesPos();
	}

	tuple<Range, string_view> getChromosomeInRefFile(string targetChromosome) {
		size_t idx = std::lower_bound(db.begin(), db.end(), ChromosomeFilePosition(std::move(targetChromosome), Range(), string_view())) - db.begin();
		assert(db[idx].name() == targetChromosome);
		return make_tuple(db[idx].range(), db[idx].slice());
	}

    string_view refFile() {
        return string_view(refFile_);
    }

private:
	class ChromosomeFilePosition {
	public:
		ChromosomeFilePosition(string&& name, Range &&range, string_view slice) : name_(name), range_(range), slice_(slice) {}

		Range range() {
			return range_;
		}

		const string& name() {
			return name_;
		}

		string_view slice() {
			return slice_;
		}

		inline bool operator<(const ChromosomeFilePosition& o) {
			assert(name_ != o.name_);
			return name_ < o.name_;
		}

	private:
		string name_;
		Range range_;		
		string_view slice_;
	};

    /**
     * given reference line (start with '>'), extract the chromosome information.
     * this is important when there is space in chromosome name. the SAM information only contain the first word.
     */
    string getChrName(string_view inputLine) {
        string name;
        for (int i = 1; i < inputLine.size(); i++)
        {
            char c = inputLine[i];
            if (isspace(c)){
                break;
            }
            name += c;
        }

        if(removedChrName) {
            if(name.find("chr") == 0) {
                name = name.substr(3);
            }
        } else if(addedChrName) {
            if(name.find("chr") != 0) {
                name = string("chr") + name;
            }
        }
        return name;
    }

	void LoadChromosomeNamesPos() {
        string_view line;
		size_t contentStart = 0; // relative offset in the file
		optional<string> chromosome = nullopt;

		// relative offsets in the file
		size_t lineStart = 0, lineEnd = 0;
		size_t lineBpCount = 0, totalBpCount = 0, location = 0;
		// totalBpCount is the total bp count of current chromosome
		// location is the start position of current chromosome, should always be 0
		// (every chromosome should map to only one record)
        for (size_t i = 0; i < refFile_.size(); ++i) {
            if (refFile_[i] == '\n') {
                line = string_view(refFile_.data() + lineStart, i - lineStart);
				if (line.front() == '>') {
					if (chromosome.has_value()) {
						db.push_back(ChromosomeFilePosition(std::move(chromosome.value()), Range(location, location + totalBpCount), string_view(refFile_.data() + contentStart, lineEnd - contentStart)));
					}
					chromosome = getChrName(line);
					contentStart = i + 1;
					totalBpCount = 0;
				} else {
					totalBpCount += lineBpCount;
				}
                lineStart = i + 1; // skip '\n'
				lineEnd = i;
				lineBpCount = 0;
            }
			if (!isspace(refFile_[i])) lineBpCount ++;
        }

        if (lineStart < refFile_.size()) {
			totalBpCount += lineBpCount;
			line = string_view(refFile_.data() + lineStart, refFile_.size() - lineStart);
			assert(line.front() != '>');
			assert(chromosome.has_value());
			db.push_back(ChromosomeFilePosition(std::move(chromosome.value()), Range(location, location + totalBpCount), string_view(refFile_.data() + contentStart, refFile_.size() - contentStart)));
        }

		sort(db.begin(), db.end());
	}

	std::vector<ChromosomeFilePosition> db;
	mio::mmap_source refFileMap;
	std::span<const char> refFile_;
    bool addedChrName = false;
    bool removedChrName = false;
};

/**
 * give a SAM line, extract the chromosome and position information.
 * return true if the SAM line is mapped. return false if SAM line is not maped.
 */
optional<tuple<string_view, long long int>> getSAMChromosomePos(string_view line) {
    string_view chr;
    long long int pos;

    int startPosition = 0;
    int endPosition = 0;
    int count = 0;

    while ((endPosition = line.find("\t", startPosition)) != string::npos) {
        if (count == 2) {
            chr = line.substr(startPosition, endPosition - startPosition);
        } else if (count == 3) {
            pos = stoll(string(line.substr(startPosition, endPosition - startPosition)));
            if (chr == "*") {
                return nullopt;
            } else {
                return make_tuple(chr, pos);
            }
        }
        startPosition = endPosition + 1;
        count++;
    }
    return nullopt;
}

string print_file_block(string_view file, string_view block) {
    auto l = block.begin() - file.begin();
    auto r = block.end() - file.begin();
    return string("[") + to_string(l) + ", " + to_string(r) + ")";
}

struct Task {
    string chromosome;
    string_view alignmentBlock;
    string_view refBlock;
};

int hisat_3n_table_2() {
    ios::sync_with_stdio(false);

	BS::thread_pool pool;

	ChromosomeDB chromosomeDB(refFileName, addedChrName, removedChrName);

	assert(!standardInMode);
	mio::mmap_source alignmentFileMap(alignmentFileName);
	string_view alignmentFile(alignmentFileMap.data(), alignmentFileMap.size());

    constexpr size_t blockLineLimit = 1048576;
    size_t blockLineCount = 0;
    vector<Task> blocks;
    string_view block = alignmentFile.substr(0, 0);
    optional<string> currentChromosome;

    // Ensure all Alignments (lines) in a block refer to the same DNA, try to split file 
    // into blocks on the DNA boundries (need the file to be sorted). 
    // Force split if block is too big that reaches blockLineLimit.
    for (auto it = LinesIter(alignmentFile); !it.empty(); ++it) {
        auto line = *it;
        if (line.front() == '@') {            
            continue;
        }
        auto coordinate = getSAMChromosomePos(line);
        if (!coordinate.has_value()) {
            continue;
        }
        auto [chromosome, location] = coordinate.value();

        if (!currentChromosome.has_value()) {
            currentChromosome = chromosome;
        }

        if (blockLineCount >= blockLineLimit || currentChromosome.value() != chromosome) {
            block = string_view(block.begin(), line.end());

            auto [range, refBlock] = chromosomeDB.getChromosomeInRefFile(string(currentChromosome.value()));
            blocks.push_back(Task { .chromosome = currentChromosome.value(), .alignmentBlock = block, .refBlock = refBlock });

            block = line.substr(0, 0);
            currentChromosome = chromosome;
        }

        blockLineCount ++;
    }

    moodycamel::ConcurrentQueue<string> outputQueue;

    // submit all blocks to thread pool
    pool.detach_sequence(0, blocks.size(), [&](size_t i) {
        const auto &b = blocks[i];
        outputQueue.enqueue(string("Worker ") + to_string(BS::this_thread::get_index().value()) + "process block " + to_string(i) + "DNA name = " + b.chromosome + ", alignmentBlock = " + print_file_block(alignmentFile, b.alignmentBlock) + ", refBlock = " + print_file_block(chromosomeDB.refFile(), b.refBlock));
    });

    while (pool.get_tasks_total() || outputQueue.size_approx()) {
        vector<string> s;
        if (outputQueue.try_dequeue_bulk(s.begin(), 10000) < 5000) {
            cout << "Sleep\n";
            this_thread::sleep_for(chrono::microseconds(10));
        }
        for (const auto &line: s) {
            cout << line << endl;
        }
    }

    return 0;
}

int main(int argc, const char** argv)
{
    int ret = 0;

    try {
        parseOptions(argc, argv);
        ret = hisat_3n_table_2();
    } catch(std::exception& e) {
        cerr << "Error: Encountered exception: '" << e.what() << "'" << endl;
        cerr << "Command: ";
        for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
        cerr << endl;
        return 1;
    } catch(int e) {
        if (e != 0) {
            cerr << "Error: Encountered internal HISAT-3N exception (#" << e << ")" << endl;
            cerr << "Command: ";
            for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
            cerr << endl;
        }
        return e;
    }

    return ret;
}
