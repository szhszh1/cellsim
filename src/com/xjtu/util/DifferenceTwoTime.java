package com.xjtu.util;

import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;

public class DifferenceTwoTime {
	/*
	 */

	public static long Difference(String t1, String t2) {
		DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		long diff = 0;
		try {
			Date d1 = df.parse(t1);
			Date d2 = df.parse(t2);

			diff = Math.abs(d1.getTime() - d2.getTime());

			
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return diff/1000;
	}
	public static long Difference1(String t1, String t2) {
		DateFormat df = new SimpleDateFormat("yyyyMMdd HH:mm:ss:SSS");
		long diff = 0;
		try {
			Date d1 = df.parse(t1);
			Date d2 = df.parse(t2);

			diff = Math.abs(d1.getTime() - d2.getTime());

			
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return diff/1000;
	}
}
