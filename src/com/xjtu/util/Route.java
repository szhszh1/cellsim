package com.xjtu.util;

import java.util.List;

import com.graphhopper.routing.AlternativeRoute;
import com.graphhopper.routing.AlternativeRoute.AlternativeInfo;
import com.graphhopper.util.GPXEntry;

public class Route {
	public List<GPXEntry> list ;
	

	public Route(List<GPXEntry> list) {
		super();
		this.list = list;
	}

	public List<GPXEntry> getList() {
		return list;
	}

	public void setList(List<GPXEntry> list) {
		this.list = list;
	}
	

	
}
