# cellsim
This is the implemetation of the code and sample data for paper "Processing Cellular Data for Similar Trajectories". The codes are implemented with MapReduce framework.

---

### Overview
1. [Requirements](#requirements)
2. [Dir](#dir)
3. [Run](#run)
4. [License](#license)

---

## 1. Requirements
The following modules are required.

* CentOS release 6.8
* Hadoop 2.6.0-cdh5.5.0
* Java 1.8 
* graphhopper-map-matching-core >= 0.8.2 (refer to https://github.com/graphhopper/map-matching)
* graph-cache (road network obtained from openstreetmap)
---

## 2. Dir
* code for map matching module

main: src\com\xjtu\mm 

mm: src\com\graphhopper\matching

* code for similarity module

main: src\com\xjtu\simi

* source code for similarity module

util: src\com\xjtu\util

* data.zip
sample data for mapmatching and similarity search

---
## 3. Run
* unzip data.zip

* run map matching
hadoop -jar mm.jar com.xjtu.mm.MM_my_hadoop_multi input output graph-cache loc_error maxnode maxroute q res

* run similar trajectory searching
hadoop -jar simi.jar com.xjtu.simi.simi_all_hadoop.java input output graph-cache maxnode globalthres query groundtruth

---
## 4. Parameter descriptions
* `input` : Path to input sequences

* `output` : Path to output results

* `graph-cache` : Path to road network (can be downloaded from openstreetmap), e.g., graph-cache

* `loc_error` : max error of the sensor e.g., 200

* `maxnode` : max visited node, e.g., 20000

* `q` : queue name in hadoop system, e.g., q1

* `res` : max java opts, e.g., 1024

* `query` : Path to query sequence, e.g., ./data/comove/2.txt

* `globalthres` : Global threshold, e.g., 300000

* `groundtruth` : Path to groundtruth sequence

---
## 4. License
The code is developed under the MPL-02.0 license.
