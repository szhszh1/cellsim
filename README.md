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

## 3. Run


---
## 4. License
The code is developed under the MPL-02.0 license.
