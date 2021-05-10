package com.xjtu.util;

public class DistanceTwoPoint {

	/*
	 
	 */
	public static double HaverSin(double theta) {
		double v = Math.sin(theta / 2);
		return v * v;
	}

	public static double EARTH_RADIUS = 6371.0;

	public static double Distance(double lat1, double lon1, double lat2,
			double lon2) {

//		lat1 = ConvertDegreesToRadians(lat1);
//		lon1 = ConvertDegreesToRadians(lon1);
//		lat2 = ConvertDegreesToRadians(lat2);
//		lon2 = ConvertDegreesToRadians(lon2);
//
//		double vLon = Math.abs(lon1 - lon2);
//		double vLat = Math.abs(lat1 - lat2);
//
//		double h = HaverSin(vLat) + Math.cos(lat1) * Math.cos(lat2)
//				* HaverSin(vLon);
//		double distance = 2 * EARTH_RADIUS * Math.asin(Math.sqrt(h));
//		return distance;
		double a, b, R;  
	    R = 6378137; // 
	    lat1 = lat1 * Math.PI / 180.0;  
	    lat2 = lat2 * Math.PI / 180.0;  
	    a = lat1 - lat2;  
	    b = (lon1 - lon2) * Math.PI / 180.0;  
	    double d;  
	    double sa2, sb2;  
	    sa2 = Math.sin(a / 2.0);  
	    sb2 = Math.sin(b / 2.0);  
	    d = 2  
	            * R  
	            * Math.asin(Math.sqrt(sa2 * sa2 + Math.cos(lat1)  
	                    * Math.cos(lat2) * sb2 * sb2));  
	    return d/1000;  
	}

	public static double ConvertDegreesToRadians(double degrees) {
		return degrees * Math.PI / 180;
	}

	public static double ConvertRadiansToDegrees(double radian) {
		return radian * 180.0 / Math.PI;
	}

}