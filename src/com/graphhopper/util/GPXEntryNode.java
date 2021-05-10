/*
 *  Licensed to GraphHopper GmbH under one or more contributor
 *  license agreements. See the NOTICE file distributed with this work for 
 *  additional information regarding copyright ownership.
 * 
 *  GraphHopper GmbH licenses this file to you under the Apache License, 
 *  Version 2.0 (the "License"); you may not use this file except in 
 *  compliance with the License. You may obtain a copy of the License at
 * 
 *       http://www.apache.org/licenses/LICENSE-2.0
 * 
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */
package com.graphhopper.util;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

import com.graphhopper.util.shapes.GHPoint;
import com.graphhopper.util.shapes.GHPoint3D;

/**
 * @author Peter Karich
 */
public class GPXEntryNode extends GHPoint3D {
    private long time;
    private int node;
    public GPXEntryNode(GHPoint p, long millis) {
        this(p.lat, p.lon, millis);
    }

    public GPXEntryNode(double lat, double lon, long millis) {
        super(lat, lon, Double.NaN);
        this.time = millis;
    }

    public GPXEntryNode(double lat, double lon, double ele, long millis) {
        super(lat, lon, ele);
        this.time = millis;
    }
    public GPXEntryNode(double lat, double lon, double ele, long millis,int node) {
        super(lat, lon, ele);
        this.time = millis;
        this.node = node;
    }
    boolean is3D() {
        return !Double.isNaN(ele);
    }

    /**
     * The time relative to the start time in milli seconds.
     */
    public long getTime() {
        return time;
    }

    public void setTime(long time) {
        this.time = time;
    }

    @Override
    public int hashCode() {
        return 59 * super.hashCode() + (int) (time ^ (time >>> 32));
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null)
            return false;

        final GPXEntryNode other = (GPXEntryNode) obj;
        return super.equals(obj);
    }
    public boolean equals1(Object obj) {
        if (obj == null)
            return false;

        final GPXEntryNode other = (GPXEntryNode) obj;
        return lat==other.lat && lon==other.lon;
    }
    public boolean equals2(Object obj) {
        if (obj == null)
            return false;

        final GPXEntryNode other = (GPXEntryNode) obj;
        return time==other.getTime() && lat==other.lat && lon==other.lon;
    }
    @Override
    public String toString() {
    	DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        return df.format(new Date(time))+","+super.lon + "," +super.lat + ","+super.ele;
    }
    public String toString1() {
    	DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        return df.format(new Date(time))+","+super.lat + "," +super.lon + ","+super.ele;
    }
}
