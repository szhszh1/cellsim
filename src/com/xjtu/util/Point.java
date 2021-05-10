package com.xjtu.util;

import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;

public class Point {
	@Override
	public String toString() {
		// TODO Auto-generated method stub
		return time+"\t"+lon+"\t"+lat;
	}
	private String time;// the time a point
	private double lon;// the coordinate of a point
	private double lat;
	
	public Point(String time, double lon, double lat) {
		super();
		this.time = time;
		this.lon = lon;
		this.lat = lat;
	}
	
	public String getTime() {
		return time;
	}
	public void setTime(String time) {
		this.time = time;
	}
	public double getLon() {
		return lon;
	}
	public void setLon(double lon) {
		this.lon = lon;
	}
	public double getLat() {
		return lat;
	}
	public void setLat(double lat) {
		this.lat = lat;
	}

	
	public boolean equals1(Point a) {
		if(a.getLat()==lat && a.getLon()==lon){
			return true;
		}else{
			return false;
		}
			
		
	}

	public boolean equals2(Point a) throws ParseException {
		DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		long t = Math.abs(df.parse(a.getTime()).getTime()-df.parse(time).getTime())/1000;
//		System.out.println(t);
		if(a.getLat()==lat && a.getLon()==lon && t<1800){
			return true;
		}else{
			return false;
		}
	}
	
	

}
