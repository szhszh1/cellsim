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
package com.xjtu.util;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

import com.graphhopper.util.shapes.GHPoint;
import com.graphhopper.util.shapes.GHPoint3D;

/**
 * @author Peter Karich
 */
public class GPXEntry {
    private long time;
    private double lat;
    private double lon;
    private double ele;

    
    
    
    
    public GPXEntry(double lat, double lon, double ele,long time) {
		super();
		this.time = time;
		this.lat = lat;
		this.lon = lon;
		this.ele = ele;
	}

	public long getTime() {
		return time;
	}

	public void setTime(long time) {
		this.time = time;
	}

	public double getLat() {
		return lat;
	}

	public void setLat(double lat) {
		this.lat = lat;
	}

	public double getLon() {
		return lon;
	}

	public void setLon(double lon) {
		this.lon = lon;
	}

	public double getEle() {
		return ele;
	}

	public void setEle(double ele) {
		this.ele = ele;
	}

	@Override
    public int hashCode() {
        return 59 * super.hashCode() + (int) (time ^ (time >>> 32));
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null)
            return false;

        final GPXEntry other = (GPXEntry) obj;
        return super.equals(obj);
    }
    public boolean equals1(Object obj) {
        if (obj == null)
            return false;

        final GPXEntry other = (GPXEntry) obj;
        return lat==other.lat && lon==other.lon;
    }
    public boolean equals2(Object obj) {
        if (obj == null)
            return false;

        final GPXEntry other = (GPXEntry) obj;
        return time==other.getTime() && lat==other.lat && lon==other.lon;
    }
    @Override
    public String toString() {
    	DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        return df.format(new Date(time))+","+lon + "," +lat + ","+ele;
    }
    public String toString1() {
    	DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        return df.format(new Date(time))+","+lat + "," +lon + ","+ele;
    }
}
