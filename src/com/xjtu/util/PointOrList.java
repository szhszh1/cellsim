package com.xjtu.util;

import java.util.List;

import com.graphhopper.util.GPXEntry;



public class PointOrList {
	public List<Route> list ;
	
	public GPXEntry p;

	
	public PointOrList(List<Route> list, GPXEntry p) {
		super();
		this.list = list;
		this.p = p;
	}

	public List<Route> getList() {
		return list;
	}

	public void setList(List<Route> list) {
		this.list = list;
	}

	public GPXEntry getP() {
		return p;
	}

	public void setP(GPXEntry p) {
		this.p = p;
	}

	@Override
	public String toString() {
		if (list == null) {
			return p.toString();
		} else {
			StringBuilder sb = new StringBuilder("");
				for (int k = 0; k < list.size(); k++) {
					for (int r = 0; r < list.get(k).getList().size(); r++) {
						sb.append(list.get(k).getList().get(r).toString());
					}
					
				}
				return sb.toString();

		}

//		return super.toString();
	}
	
	
	
	
	
}
