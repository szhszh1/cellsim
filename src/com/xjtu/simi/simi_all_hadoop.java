package com.xjtu.simi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.filecache.DistributedCache;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.Reducer.Context;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.MultipleOutputs;

import com.graphhopper.GraphHopper;
import com.graphhopper.matching.MapMatching_my;
import com.graphhopper.reader.osm.GraphHopperOSM;
import com.graphhopper.routing.AlgorithmOptions;
import com.graphhopper.routing.Dijkstra;
import com.graphhopper.routing.util.DefaultEdgeFilter;
import com.graphhopper.routing.util.EdgeFilter;
import com.graphhopper.routing.util.FlagEncoder;
import com.graphhopper.routing.util.HintsMap;
import com.graphhopper.routing.util.TraversalMode;
import com.graphhopper.routing.weighting.ShortestWeighting;
import com.graphhopper.routing.weighting.Weighting;
import com.graphhopper.storage.GraphHopperStorage;
import com.graphhopper.storage.index.LocationIndex;
import com.graphhopper.storage.index.QueryResult;
import com.graphhopper.util.CmdArgs;
import com.graphhopper.util.GPXEntry;
import com.graphhopper.util.Parameters;
import com.xjtu.util.PointOrList;
import com.xjtu.util.Route;
import com.xjtu.util.DistanceTwoPoint;

/**
 * 
 *  
 *  1. my simi
 *  2. vldb simi
 *  
 *  
 *  map 
 *  	simi_my       user,simi
 *  	simi_vldb	  user,simi
 *  
 *  reduce
 *  	user 1,5,10
 *  
 *  @param
 *  	input			/
 *  	output			/
 *  	graph			graph-cache
 *      gpsAccuracy     200
 * 	 	maxNode			20000
 *		global			300000
 *		query           query traj path
 *		groundtruth     groundtruth path
 *  	
 */
public class simi_all_hadoop {
	public static class SimiMap extends Mapper<LongWritable, Text, Text, Text> {
		GraphHopper hopper;
		DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		Map<String,String> map = new HashMap<String,String>();
		protected void setup(Context context) throws IOException,
				InterruptedException {
			Configuration conf = context.getConfiguration();
			// query trajectory  
			String query = conf.get("query");

			FileSystem data = FileSystem.get(conf);

			FSDataInputStream in1 = data.open(new Path(query));
			BufferedReader br1 = new BufferedReader(new InputStreamReader(in1));
			while (true) {
				String line = null;
				line = br1.readLine();
				if (line == null) {
					break;
				}
				map.put("query", line);
			}
			
			
			//2
			Path[] localCacheFiles = DistributedCache
					.getLocalCacheFiles(context.getConfiguration());
			String s = "action=import datasource=./some-dir/osm-file.pbf vehicle=car";
			String[] arg = s.split(" ");
			CmdArgs args = CmdArgs.read(arg);
			args.put("action", "import");
			args.put("graph.location", localCacheFiles[0].toUri().getPath()
					.toString());
			args.put("vehicles", "car");
			args.put("graph.flag_encoders", "car");
			args.put("datareader.file", localCacheFiles[0].toUri().getPath()
					.toString());
			args.put("prepare.min_one_way_network_size", 200);
			hopper = new GraphHopperOSM().init(args);
			hopper.getCHFactoryDecorator().setEnabled(false);
			hopper.importOrLoad();
		}

		public void map(LongWritable key, Text value, Context context)
				throws IOException, InterruptedException {
//			try{
				Configuration conf = context.getConfiguration();
				long timeThres = conf.getLong("timeThres", 300000);
				double spaceThres = conf.getDouble("spaceThres", 0.1);
				double geoThres = conf.getDouble("geoThres", 0.05);
				double gpsAccuracy = conf.getDouble("gpsAccuracy", 200);
				int maxroute = conf.getInt("maxroute", 10);
				int maxNode = conf.getInt("maxNode", 5000);
				
				FlagEncoder firstEncoder = hopper.getEncodingManager()
						.fetchEdgeEncoders().get(0);
				
				AlgorithmOptions opts = AlgorithmOptions
						.start()
						.algorithm(Parameters.Algorithms.DIJKSTRA)
						.traversalMode(TraversalMode.EDGE_BASED_2DIR)
						.weighting(new ShortestWeighting(firstEncoder))
						.maxVisitedNodes(maxNode)
						.hints(new HintsMap().put("weighting", "shortest").put(
								"vehicle", firstEncoder.toString())).build();
				
				MapMatching_my mapMatching_my = new MapMatching_my(hopper, opts);
				mapMatching_my.setTransitionProbabilityBeta(0.00959442);
				mapMatching_my.setMeasurementErrorSigma(gpsAccuracy);
				GraphHopperStorage ghs = hopper.getGraphHopperStorage();
				Weighting weighting = opts.getWeighting();
				TraversalMode traversalMode = opts.getTraversalMode();
				EdgeFilter edgeFilter = new DefaultEdgeFilter(opts.getWeighting().getFlagEncoder());
				
				String query = map.get("query");
				String[] refer_user_traj = query.split("\t",-1);
				String[] wait_user_traj = value.toString().split("\t",-1);
				String referuser = refer_user_traj[0];
				String waituser = wait_user_traj[0];
				List<PointOrList> referTraj = tranLineToPointOrList(refer_user_traj[1]);
				List<PointOrList> waitTraj = tranLineToPointOrList(wait_user_traj[1]);
				
				if(referTraj.size()>0 && waitTraj.size()>0){
					boolean bool = global(referTraj, waitTraj, timeThres);
					if (bool) {
						long start = getStartTime(referTraj);
						long end = getEndTime(referTraj);

						List<Route> referAll = transToMultiRoute(maxroute,
								referTraj);
						List<Route> waitAll1 = transToMultiRoute(maxroute,
								waitTraj);
						List<Route> waitAll = getSubWait(waitAll1, start, end);
						if(waitAll.size()>0){
							List<GPXEntry> listrefer = getObserGPXEntry(referTraj);
							List<GPXEntry> listwait = getObserGPXEntry(waitTraj);
							
							double simi_my = simiFunc_my_all(referAll, waitAll,
									listrefer, timeThres, spaceThres, geoThres);
							double simi_vldb = simiFunc_vldb_all(referAll, waitAll, ghs,
									weighting, traversalMode, hopper, edgeFilter);
							if(simi_my>0){
								context.write(new Text(referuser+",simi_my"), new Text(waituser+","+String.valueOf(simi_my)));
							}
							if(simi_vldb>0){
								context.write(new Text(referuser+",simi_vldb"), new Text(waituser+","+String.valueOf(simi_vldb)));
							}
						}
						
					}
				}
				
//			}catch (Exception e){
//				
//			}
			
		}
	}
	
	
	public static class SimiReduce extends Reducer<Text, Text, Text, Text> {
		private MultipleOutputs output;

        @Override
        protected void setup(Context context
        ) throws IOException, InterruptedException {
            output = new MultipleOutputs(context);
        }
        
        public void reduce(Text key, Iterable<Text> values, Context context)
				throws IOException, InterruptedException {
			Map<String, Double> map = new TreeMap<String, Double>();
			for (Text i : values) {
				String[] seg = i.toString().split(",");
				map.put(seg[0], Double.parseDouble(seg[1]));
			}
			map = sortMap(map);

			Iterator<Entry<String, Double>> it = map.entrySet().iterator();
			while (it.hasNext()) {
				Entry<String, Double> entry = it.next();
				String user = entry.getKey();
				if(key.toString().contains("simi_my")){
					output.write(new Text(user),
							new Text(String.valueOf(entry.getValue())), "simi_my");
				}
				
				if(key.toString().contains("simi_vldb")){
					output.write(new Text(user),
							new Text(String.valueOf(entry.getValue())), "simi_vldb");
				}
				
			}

		}
        @Override
        protected void cleanup(Context context
        ) throws IOException, InterruptedException {
            output.close();
        }
	}

	
	
	
	public static void process(String query,String company,Map<String, Map<String, String>> map) {
		if (map.containsKey(query)) {
			Map<String, String> m = map.get(query);
			m.put(company, "");
			map.put(query, m);
		} else {
			Map<String, String> m = new HashMap<String, String>();
			m.put(company, "");
			map.put(query, m);
		}
	}
	
	public static Map<String,Double> sortMap(Map<String,Double> oldMap) {
		ArrayList<Map.Entry<String,Double>> list = new ArrayList<Map.Entry<String,Double>>(
				oldMap.entrySet());
		Collections.sort(list, new Comparator<Map.Entry<String,Double>>() {

			@Override
			public int compare(Map.Entry<String,Double> o1, Map.Entry<String,Double> o2) {  
                return o2.getValue().compareTo(o1.getValue());  
           }  
		});
		Map<String,Double> newMap = new LinkedHashMap<String,Double>();
		for (int i = 0; i < list.size(); i++) {
			newMap.put(list.get(i).getKey(), list.get(i).getValue());
		}
		return newMap;
	}
	private static double simiFunc_vldb_all(List<Route> referAll, List<Route> waitAll,
			GraphHopperStorage ghs, Weighting weighting,
			TraversalMode traversalMode, GraphHopper hopper,
			EdgeFilter edgeFilter) {
		double simiGlobal = 0;
		for (int index_refer = 0; index_refer < referAll.size(); index_refer++) {
			double simiLocal = 0;
			List<GPXEntry> referListCurrent = referAll.get(index_refer)
					.getList();
			for (int index_wait = 0; index_wait < waitAll.size(); index_wait++) {
				List<GPXEntry> waitListCurrent = waitAll.get(index_wait)
						.getList();
				List<GPXEntry> listAll = new ArrayList<GPXEntry>();
				listAll.addAll(waitListCurrent);
				listAll.addAll(referListCurrent);
				double simi_cur = simiFunc_vldb(referListCurrent,
						waitListCurrent, ghs, weighting, traversalMode, hopper,
						edgeFilter);
				if (simi_cur > simiLocal) {
					simiLocal = simi_cur;
				}
			}
			if (simiLocal > simiGlobal) {
				simiGlobal = simiLocal;
			}
		}
		return simiGlobal/2;
	}

	private static double simiFunc_vldb(List<GPXEntry> referTraj1,
			List<GPXEntry> waitTraj1, GraphHopperStorage ghs,
			Weighting weighting, TraversalMode traversalMode,
			GraphHopper hopper, EdgeFilter edgeFilter) {
		
		List<GPXEntry> referTraj = transTo1Traj(referTraj1);
		List<GPXEntry> waitTraj = transTo1Traj(waitTraj1);
		double timeSimi1 = calTimeSimi(referTraj, waitTraj);
		double timeSimi2 = calTimeSimi(waitTraj, referTraj);
		double spaceSimi1 = calSpaceSimi(referTraj, waitTraj, ghs, weighting,
				traversalMode, hopper, edgeFilter);
		
		double spaceSimi2 = calSpaceSimi(waitTraj, referTraj, ghs, weighting,
				traversalMode, hopper, edgeFilter);
		

		double simi = 0.5 * ((timeSimi1) + (timeSimi2)) + 0.5
				* ((spaceSimi1) + (spaceSimi2));
		return simi;
	}
	private static List<GPXEntry> transTo1Traj(List<GPXEntry> traj) {
		List<GPXEntry> list = new ArrayList<GPXEntry>();
		for(int i=0;i<traj.size();i++){
			if(traj.get(i).getEle()==1){
				list.add(traj.get(i));
			}
		}
		return list;
	}

	private static double calSpaceSimi(List<GPXEntry> referTraj,
			List<GPXEntry> waitTraj, GraphHopperStorage ghs,
			Weighting weighting, TraversalMode traversalMode,
			GraphHopper hopper, EdgeFilter edgeFilter) {
		double spaceSimi = 0;
		for (int i = 0; i < referTraj.size(); i++) {
			double d = minDis(referTraj.get(i), waitTraj, ghs, weighting,
					traversalMode, hopper, edgeFilter);

			spaceSimi = spaceSimi + Math.exp(-d / 1000);
			
		}
		return spaceSimi / referTraj.size();
	}

	private static double minDis(GPXEntry g, List<GPXEntry> waitTraj,
			GraphHopperStorage ghs, Weighting weighting,
			TraversalMode traversalMode, GraphHopper hopper,
			EdgeFilter edgeFilter) {
		LocationIndex loc = hopper.getLocationIndex();

		double min = Double.MAX_VALUE;
		for (int i = 0; i < waitTraj.size(); i++) {

			GPXEntry g1 = waitTraj.get(i);
			Dijkstra dij = new Dijkstra(ghs, weighting, traversalMode);
			double d = Double.MAX_VALUE;
			QueryResult from = loc.findClosest(g.getLon(), g.getLat(), edgeFilter);
			QueryResult to = loc.findClosest(g1.getLon(), g1.getLat(), edgeFilter);
			
			try {
				d = dij.calcPath(from.getClosestNode(), to.getClosestNode())
						.getDistance();
			} catch (Exception e) {
			}
			if (d <= min) {
				min = d;
			}
		}
		return min;
	}

	private static double calTimeSimi(List<GPXEntry> referTraj,
			List<GPXEntry> waitTraj) {
		double timeSimi = 0;
		for (int i = 0; i < referTraj.size(); i++) {
			GPXEntry g = referTraj.get(i);
			int d = minTime(g, waitTraj);
			timeSimi = timeSimi + Math.exp(-d / 60000);
		}
		return timeSimi / referTraj.size();
	}

	private static int minTime(GPXEntry g, List<GPXEntry> waitTraj) {
		int min = Integer.MAX_VALUE;
		for (int i = 0; i < waitTraj.size(); i++) {
			GPXEntry g1 = waitTraj.get(i);
			int minus = (int) Math.abs(g.getTime() - g1.getTime());
			if (minus < min) {
				min = minus;
			}
		}
		return min;
	}

	private static double simiFunc_my_all(List<Route> referAll,
			List<Route> waitAll, List<GPXEntry> listrefer, long timeThres,
			double spaceThres,double geoThres) {
		double simiGlobal = 0;
		int c = 0;
		for (int index_refer = 0; index_refer < referAll.size(); index_refer++) {

			double simiLocal = 0;
			List<GPXEntry> referListCurrent = referAll.get(index_refer)
					.getList();
			for (int index_wait = 0; index_wait < waitAll.size(); index_wait++) {
				List<GPXEntry> waitListCurrent = waitAll.get(index_wait)
						.getList();
				boolean bool1 = Time(referListCurrent, waitListCurrent,
						timeThres);

				boolean bool2 = Space(referListCurrent, waitListCurrent,
						spaceThres);

				if (bool1 && bool2) {
					double simi_cur = simiFunc_my(listrefer, referListCurrent,
							waitListCurrent, geoThres);
					// System.out.println(simi_cur);
					// System.out.println(referAll.get(c).getList().size());

					if (simi_cur > simiLocal) {
						simiLocal = simi_cur;
					}
				}
				if (simiLocal > simiGlobal) {
					simiGlobal = simiLocal;
					c = index_refer;
				}
			}
		}

		return simiGlobal / (referAll.get(c).getList().size() - 1);
	}
	private static boolean Space(List<GPXEntry> referListCurrent,
			List<GPXEntry> waitListCurrent, double spaceThres) {
		GPXEntry g1r = referListCurrent.get(0);
		GPXEntry g2r = referListCurrent.get(referListCurrent.size() - 1);
		GPXEntry g1w = waitListCurrent.get(0);
		GPXEntry g2w = waitListCurrent.get(waitListCurrent.size() - 1);

		double d1 = DistanceTwoPoint.Distance(g1r.getLat(), g1r.getLon(),
				g1w.getLat(), g1w.getLon());
		double d2 = DistanceTwoPoint.Distance(g2r.getLat(), g2r.getLon(),
				g2w.getLat(), g2w.getLon());
		if (d1 < spaceThres && d2 < spaceThres) {
			return true;
		}
		return false;

	}

	private static boolean Time(List<GPXEntry> referListCurrent,
			List<GPXEntry> waitListCurrent, long timeThres) {
		GPXEntry g1r = referListCurrent.get(0);
		GPXEntry g2r = referListCurrent.get(referListCurrent.size() - 1);
		GPXEntry g1w = waitListCurrent.get(0);
		GPXEntry g2w = waitListCurrent.get(waitListCurrent.size() - 1);

		if (Math.abs(g1r.getTime() - (g1w.getTime())) < timeThres
				&& Math.abs(g2r.getTime() - (g2w.getTime())) < timeThres) {
			return true;
		}
		return false;

	}
	private static double simiFunc_my(List<GPXEntry> listrefer,
			List<GPXEntry> referListCurrent, List<GPXEntry> waitListCurrent,double geoThres) {
		double simi = 0;
		for (int i = 0; i < listrefer.size() - 1; i++) {
			GPXEntry g1 = listrefer.get(i);
			GPXEntry g2 = listrefer.get(i + 1);

			List<GPXEntry> referList = subListByPoint(referListCurrent, g1, g2);
			int loczuo = minTime_zuo(g1.getTime(), waitListCurrent);
			int locyou = minTime_you(g2.getTime(), waitListCurrent);
			List<GPXEntry> waitList = subListByIndex(waitListCurrent, loczuo,
					locyou);
			//
			simi += union(referList, waitList,geoThres);
			// simi += union(referListCurrent, waitListCurrent);
		}
		return simi;
	}

	private static double union(List<GPXEntry> referList,
			List<GPXEntry> waitList,double geoThres) {
		List<GPXEntry> realA = new ArrayList<GPXEntry>(referList);
		List<GPXEntry> realB = new ArrayList<GPXEntry>(waitList);
		Set<String> setA = new HashSet<String>();
		Set<String> setB = new HashSet<String>();
		for (int i = 0; i < realA.size() - 1; i++) {
			setA.add(realA.get(i).lat + "," + realA.get(i).lon + "|"
					+ realA.get(i + 1).lat + "," + realA.get(i + 1).lon);
		}
		for (int i = 0; i < realB.size() - 1; i++) {
			setB.add(realB.get(i).lat + "," + realB.get(i).lon + "|"
					+ realB.get(i + 1).lat + "," + realB.get(i + 1).lon);
		}
		int s = geoNear(setA, setB,geoThres);
		return s;
		// setA.retainAll(setB);
		// return setA.size();
	}

	private static int geoNear(Set<String> setA, Set<String> setB,double geoThres) {
		int size = 0;
		Iterator<String> it = setA.iterator();
		while (it.hasNext()) {
			String str = it.next();
			Iterator<String> it1 = setB.iterator();
			while (it1.hasNext()) {
				String str1 = it1.next();
				double d = geo(str, str1);
				if (d < 0.005) {
					size++;
				}

			}
		}
		if (size > setA.size()) {
			size = setA.size();
		}
		return size;
	}

	private static double geo(String str, String str1) {
		String[] x = str.split("\\|");
		String[] y = str1.split("\\|");
		String[] ll11 = x[0].split(",");
		String[] ll12 = x[1].split(",");
		String[] ll21 = y[0].split(",");
		String[] ll22 = y[1].split(",");

		double d = DistanceTwoPoint.Distance(Double.parseDouble(ll11[0]),
				Double.parseDouble(ll11[1]), Double.parseDouble(ll21[0]),
				Double.parseDouble(ll21[1]));
		double d1 = DistanceTwoPoint.Distance(Double.parseDouble(ll12[0]),
				Double.parseDouble(ll12[1]), Double.parseDouble(ll22[0]),
				Double.parseDouble(ll22[1]));
		return Math.min(d, d1);
	}
	private static List<Route> getSubWait(List<Route> traj, long start, long end) {
		List<Route> result = new ArrayList<Route>();
		for (int i = 0; i < traj.size(); i++) {
			Route route = sub(traj.get(i), start, end);
			if(route.getList().size()>0){
				result.add(route);
			}
			
		}
		return result;
	}

	private static Route sub(Route route, long start, long end) {
		List<GPXEntry> list = new ArrayList<GPXEntry>();
		int start_index = minTime_zuo(start, route.getList());
		int end_index = minTime_you(end, route.getList());
		if(start_index<end_index){
			for (int i = start_index; i <= end_index; i++) {
				list.add(route.getList().get(i));
			}
		}
		return new Route(list);
	}

	private static long getEndTime(List<PointOrList> traj) {
		int len = traj.size() - 1;
		if (traj.get(len).getList() != null) {
			List<Route> l = traj.get(len).getList();
			return l.get(len).getList().get(len).getTime();
		} else {
			return traj.get(len).getP().getTime();
		}
	}

	private static long getStartTime(List<PointOrList> traj) {
		if (traj.get(0).getList() != null) {
			List<Route> l = traj.get(0).getList();
			return l.get(0).getList().get(0).getTime();
		} else {
			return traj.get(0).getP().getTime();
		}
	}

	private static boolean global(List<PointOrList> referTraj,
			List<PointOrList> waitTraj, long timeThres) {
		GPXEntry g1r = getPointOrListElement1(referTraj.get(0));
		GPXEntry g2r = getPointOrListElement(referTraj
				.get(referTraj.size() - 1));
		GPXEntry g1w = getPointOrListElement1(waitTraj.get(0));
		GPXEntry g2w = getPointOrListElement(waitTraj.get(waitTraj.size() - 1));

		if (g1r.getTime() >= (g1w.getTime() - timeThres)
				&& g2r.getTime() <= (g2w.getTime() + timeThres)) {
			return true;
		}
		return false;
	}

	private static GPXEntry getPointOrListElement(PointOrList pointOrList) {
		if (pointOrList.getList() != null) {
			List<GPXEntry> l = pointOrList.getList().get(0).getList();
			GPXEntry g = l.get(l.size() - 1);
			return g;
		} else {
			GPXEntry g = pointOrList.getP();
			return g;
		}
	}

	private static GPXEntry getPointOrListElement1(PointOrList pointOrList) {
		if (pointOrList.getList() != null) {
			List<GPXEntry> l = pointOrList.getList().get(0).getList();
			GPXEntry g = l.get(0);
			return g;
		} else {
			GPXEntry g = pointOrList.getP();
			return g;
		}
	}

	private static List<GPXEntry> subListByIndex(List<GPXEntry> list,
			int loczuo, int locyou) {
		List<GPXEntry> result = new ArrayList<GPXEntry>();
		for (int i = loczuo; i <= locyou; i++) {
			result.add(list.get(i));
		}
		return result;
	}

	private static List<GPXEntry> subListByPoint(List<GPXEntry> list,
			GPXEntry g1, GPXEntry g2) {
		List<GPXEntry> result = new ArrayList<>();
		boolean b = false;
		for (int i = 0; i < list.size(); i++) {
			if ((list.get(i).lat == g1.lat && list.get(i).lon == g1.lon && list
					.get(i).getTime() == g1.getTime()) || b) {
				b = true;
				result.add(list.get(i));
			}
			if (list.get(i).lat == g2.lat && list.get(i).lon == g2.lon
					&& list.get(i).getTime() == g2.getTime()) {
				b = false;
			}
		}
		return result;
	}

	private static double union(List<GPXEntry> referList,
			List<GPXEntry> waitList) {
		List<GPXEntry> realA = new ArrayList<GPXEntry>(referList);
		List<GPXEntry> realB = new ArrayList<GPXEntry>(waitList);
		Set<String> setA = new HashSet<String>();
		Set<String> setB = new HashSet<String>();
		for (int i = 0; i < realA.size() - 1; i++) {
			setA.add(realA.get(i).lat + ","+realA.get(i).lon
					+ realA.get(i+1).lat + ","+realA.get(i+1).lon);
		}
		for (int i = 0; i < realB.size() - 1; i++) {
			setB.add(realB.get(i).lat + ","+realB.get(i).lon
					+ realB.get(i+1).lat + ","+realB.get(i+1).lon);
		}
		setA.retainAll(setB);
		return setA.size();
	}




	private static List<Route> transToMultiRoute(int count, List<PointOrList> l) {
		List<Route> ll = new ArrayList<Route>(count);
		for (int i = 0; i < count; i++) {
			List<GPXEntry> list = new ArrayList<GPXEntry>();
			ll.add(new Route(list));
		}

		for (int i = 0; i < l.size(); i++) {
			if (l.get(i).getList() != null) {

				int size = l.get(i).getList().size();

				int k = 0;
				for (int j = 0; j < ll.size(); j++) {

					if (k < size) {
						List<GPXEntry> gg = l.get(i).getList().get(k).getList();
						for (int r = 0; r < gg.size(); r++) {

							ll.get(j).getList().add(gg.get(r));

						}
						k++;
					} else {
						k = 0;
					}

				}

			} else {

				for (int j = 0; j < ll.size(); j++) {
					ll.get(j).getList().add(l.get(i).getP());
				}
			}

		}
		return ll;

	}

	private static int minTime(long time, List<GPXEntry> listwait) {
		int min = Integer.MAX_VALUE;
		int c = 0;
		for (int i = 0; i < listwait.size(); i++) {
			GPXEntry g1 = listwait.get(i);
			int minus = (int) Math.abs(time - g1.getTime());
			if (minus < min) {
				min = minus;
				c = i;
			}
		}
		return c;

	}
	private static int minTime_zuo(long time, List<GPXEntry> listwait) {
		int c = 0;
		List<GPXEntry> list = new ArrayList<GPXEntry>();
		for (int i = 0; i < listwait.size(); i++) {
			GPXEntry g1 = listwait.get(i);
			int minus = (int) Math.abs(time - g1.getTime());
			
			if(minus<120000){
				list.add(g1);
			}
		}
		if(list.size()>0){
			GPXEntry gx = list.get(0);
			for (int i = 0; i < listwait.size(); i++) {
				if(listwait.get(i).getTime()==gx.getTime() && listwait.get(i).lat==gx.lat && listwait.get(i).lon==gx.lon){
					return i;
				}
			}
		}else{
			c=minTime(time, listwait);
			return c;
		}
		return c;
	}
	private static int minTime_you(long time, List<GPXEntry> listwait) {
		int c = 0;
		List<GPXEntry> list = new ArrayList<GPXEntry>();
		for (int i = 0; i < listwait.size(); i++) {
			GPXEntry g1 = listwait.get(i);
			int minus = (int) Math.abs(time - g1.getTime());
			if(minus<300000){
				list.add(g1);
			}
		}
//		System.out.println(list.size());
		if(list.size()>0){
			GPXEntry gx = list.get(list.size()-1);
			
			for (int i = 0; i < listwait.size(); i++) {
				
				if(listwait.get(i).getTime()==gx.getTime() && listwait.get(i).lat==gx.lat && listwait.get(i).lon==gx.lon){
					
					return i;
				}
			}
		}else{
			c=minTime(time, listwait);
			return c;
		}
		return c;
	}
	private static List<GPXEntry> getObserGPXEntry(List<PointOrList> referTraj) {
		List<GPXEntry> listrefer = new ArrayList<GPXEntry>();
		Map<Long, String> map = new HashMap<Long, String>();
		for (int i = 0; i < referTraj.size(); i++) {
			PointOrList pl = referTraj.get(i);
			if (pl.getList() == null) {
				if (pl.getP().getEle() == 1) {
					if (!map.containsKey(pl.getP().getTime())) {
						map.put(pl.getP().getTime(), "");
						listrefer.add(pl.getP());
					}
				}
			} else {
				List<GPXEntry> l = pl.getList().get(0).getList();
				if (l.get(0).getEle() == 1) {
					if (!map.containsKey(l.get(0).getTime())) {
						map.put(l.get(0).getTime(), "");
						listrefer.add(l.get(0));
					}
				}
				if (l.get(l.size() - 1).getEle() == 1) {
					if (!map.containsKey(l.get(l.size() - 1).getTime())) {
						map.put(l.get(l.size() - 1).getTime(), "");
						listrefer.add(l.get(l.size() - 1));
					}
				}
			}
		}
		return listrefer;
	}

	/**
	 * 
	 * @param traj
	 * @return List<PointOrList>
	 * 
	 *     
	 * 
	 * 
	 */


	private static int getNumCombination(List<PointOrList> traj) {
		int count = 1;
		for (int i = 0; i < traj.size(); i++) {
			if (traj.get(i).getList() != null) {
				count = count * traj.get(i).getList().size();
			}
		}
		return count;
	}
	/**
	 * 
	 * @param traj
	 * @return List<PointOrList>
	 * 
	 * 
	 * 
	 */

	private static List<PointOrList> tranLineToPointOrList(String traj) {

		DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		List<PointOrList> refer = new ArrayList<PointOrList>();
		if (!(traj.indexOf("&") == -1)) {
			String[] pl = traj.split("&", -1);
			for (int i = 0; i < pl.length; i++) {

				if (!(pl[i].indexOf("|") == -1) || !(pl[i].indexOf("@") == -1)) {

					if (!(pl[i].indexOf("|") == -1)) {

						String[] alt = pl[i].split("\\|", -1);

						GPXEntry gx = new GPXEntry(0, 0, 0, 0);
						List<Route> route = new ArrayList<Route>();
						int len = alt.length;
						if (len >= 5) {
							len = 5;
						}
						for (int j = 0; j < len; j++) {

							if (!(alt[j].indexOf("@") == -1)) {
								String[] location = alt[j].split("@", -1);
								List<GPXEntry> list = new ArrayList<GPXEntry>();
								for (int k = 0; k < location.length; k++) {
									String[] seg = location[k].split(",");
									double lat = Double.parseDouble(seg[2]);
									double lon = Double.parseDouble(seg[1]);

									double ele = Double.parseDouble(seg[3]);
									long time = 0;
									try {
										time = df.parse(seg[0]).getTime();
									} catch (ParseException e) {
									}
									list.add(new GPXEntry(lat, lon, ele, time));
								}
								Route r = new Route(list);
								route.add(r);

							}
						}
						PointOrList p = new PointOrList(route, gx);
						refer.add(p);

					} else {
						GPXEntry gx = new GPXEntry(0, 0, 0, 0);
						List<Route> route = new ArrayList<Route>();
						if (!(pl[i].indexOf("@") == -1)) {
							String[] location = pl[i].split("@", -1);
							List<GPXEntry> list = new ArrayList<GPXEntry>();
							for (int k = 0; k < location.length; k++) {

								String[] seg = location[k].split(",");
								double lat = Double.parseDouble(seg[2]);
								double lon = Double.parseDouble(seg[1]);
								double ele = Double.parseDouble(seg[3]);
								long time = 0;
								try {
									time = df.parse(seg[0]).getTime();
								} catch (ParseException e) {
								}
								list.add(new GPXEntry(lat, lon, ele, time));
							}
							Route r = new Route(list);
							route.add(r);
						}
						PointOrList p = new PointOrList(route, gx);
						refer.add(p);
					}

				} else {
					List<Route> route = null;
					String[] seg = pl[i].split(",");
					double lat = Double.parseDouble(seg[2]);
					double lon = Double.parseDouble(seg[1]);
					double ele = Double.parseDouble(seg[3]);
					long time = 0;
					try {
						time = df.parse(seg[0]).getTime();
					} catch (ParseException e) {
					}
					GPXEntry gx = new GPXEntry(lat, lon, ele, time);
					PointOrList p = new PointOrList(route, gx);
					refer.add(p);
				}
			}

		}
		return refer;
	}
	
	public static void main(String[] args) throws Exception {
		Configuration conf = new Configuration();
		conf.set("mapred.textoutputformat.ignoreseparator", "true");
		conf.set("mapred.textoutputformat.separator", "\t");
		conf.set("mapreduce.job.priority", "VERY_HIGH");
		conf.set("mapred.child.java.opts", "-Xmx2048m");
		conf.set("mapreduce.map.failures.maxpercent", "1");
		conf.set("mapred.task.timeout", "0");
		conf.set("gpsAccuracy", args[3]);
		conf.set("maxNode", args[4]);
		// 300000
		conf.set("timeThres", args[5]);
		conf.set("spaceThres", args[6]);
		conf.set("geoThres", args[7]);
		// query traj
		conf.set("query", args[8]);
		conf.set("maxroute", args[9]);
		//groundtruth
		// query \t wait
//		conf.set("groundtruth", args[7]);
		
		
		Path vocabularyPath = new Path(args[2]);
		DistributedCache.addCacheFile(vocabularyPath.toUri(), conf);
		Job job = new Job(conf, "Simi");
		job.setJarByClass(simi_all_hadoop.class);
		job.setMapperClass(SimiMap.class);
		job.setReducerClass(SimiReduce.class);
		job.setNumReduceTasks(2);
		job.setOutputKeyClass(Text.class);
		job.setOutputValueClass(Text.class);
		FileInputFormat.addInputPath(job, new Path(args[0]));
		FileOutputFormat.setOutputPath(job, new Path(args[1]));
		System.exit(job.waitForCompletion(true) ? 0 : 1);
	}

}
