package com.xjtu.mm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URI;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.filecache.DistributedCache;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;

import com.graphhopper.GraphHopper;
import com.graphhopper.matching.EdgeMatch;
import com.graphhopper.matching.MapMatching_my;
import com.graphhopper.matching.MatchResult;
import com.graphhopper.reader.osm.GraphHopperOSM;
import com.graphhopper.routing.AlgorithmOptions;
import com.graphhopper.routing.AlternativeRoute;
import com.graphhopper.routing.Dijkstra;
import com.graphhopper.routing.AlternativeRoute.AlternativeInfo;
import com.graphhopper.routing.util.FlagEncoder;
import com.graphhopper.routing.util.HintsMap;
import com.graphhopper.routing.util.TraversalMode;
import com.graphhopper.routing.weighting.ShortestWeighting;
import com.graphhopper.routing.weighting.Weighting;
import com.graphhopper.storage.GraphHopperStorage;
import com.graphhopper.util.CmdArgs;
import com.graphhopper.util.GPXEntry;
import com.graphhopper.util.Parameters;
import com.graphhopper.util.PointList;
import com.graphhopper.util.shapes.GHPoint3D;
import com.xjtu.util.DistanceTwoPoint;
import com.xjtu.util.Point;
import com.xjtu.util.PointOrList;
import com.xjtu.util.Route;

/**
 * 
 * @param input ***
 * 			output ***
 * 			 graph graph-cache 
 * 			gpsAccuracy 200 
 * 			column 3
 * 			maxNode  20000
 * 
 */
public class MM_my_hadoop_multi {
	public static class MMMap extends Mapper<LongWritable, Text, Text, Text> {
		GraphHopper hopper;
		DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		protected void setup(Context context) throws IOException,
				InterruptedException {
			Path[] localCacheFiles = DistributedCache
					.getLocalCacheFiles(context.getConfiguration());
			// URI[] localCacheFiles = context.getCacheFiles();
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
			Text word = new Text();
			Text val = new Text();
			try {
				Configuration conf = context.getConfiguration();
				double gpsAccuracy = conf.getDouble("gpsAccuracy", 200);
				int maxNode = conf.getInt("maxNode", 20000);
				int maxPath = conf.getInt("maxpath", 10);
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

				// mapmatching 
				MapMatching_my mapMatching_tits = new MapMatching_my(hopper, opts);
				mapMatching_tits.setTransitionProbabilityBeta(0.0959442);
				mapMatching_tits.setMeasurementErrorSigma(gpsAccuracy);

				GraphHopperStorage ghs = hopper.getGraphHopperStorage();
				Weighting weighting = opts.getWeighting();
				TraversalMode traversalMode = opts.getTraversalMode();

				/**
				 * Part 3 mm
				 */
				String[] sg = value.toString().split("\t");
				DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
				List<PointOrList> result = mm_cellular_my(sg[1],
						mapMatching_tits, ghs, weighting, traversalMode,maxPath);
				StringBuilder result1 = new StringBuilder(sg[0]);
				result1.append("\t");
				for (int i = 0; i < result.size(); i++) {
					PointOrList pl = result.get(i);
					if (pl.getList() == null) {

						result1.append(df.format(new Date(pl.getP().getTime())))
								.append(",");
						result1.append(pl.getP().getLat()).append(",");
						result1.append(pl.getP().getLon()).append(",");
						result1.append(pl.getP().getEle());
						result1.append("&");
					} else {
						List<Route> list1 = pl.getList();
						for (int k = 0; k < list1.size(); k++) {
							for (int r = 0; r < list1.get(k).getList().size(); r++) {
								result1.append(
										df.format(new Date(list1.get(k)
												.getList().get(r).getTime())))
										.append(",");
								result1.append(
										list1.get(k).getList().get(r).getLat())
										.append(",");
								result1.append(
										list1.get(k).getList().get(r).getLon())
										.append(",");
								result1.append(
										list1.get(k).getList().get(r).getEle())
										.append("@");
							}
							result1.deleteCharAt(result1.length() - 1);
							result1.append("|");
						}
						result1.deleteCharAt(result1.length() - 1);
						result1.append("&");
					}

				}
				word.set(result1.deleteCharAt(result1.length() - 1).toString());
				val.set("");
				context.write(word, val);
			} catch (Exception e) {
			}
		}
	}

	private static List<PointOrList> mm_cellular_my(String cellPath,
			MapMatching_my mapMatching_my, GraphHopperStorage ghs,
			Weighting weighting, TraversalMode traversalMode, int maxPath)
			throws ParseException, NumberFormatException, IOException {

		List<Point> list_my = readLine(cellPath);

		List<GPXEntry> list_pre_my_gpx1 = null;
		try {
			list_pre_my_gpx1 = PointtoGPXEntry(list_my);
		} catch (Exception e) {
		}
		Map<Long, String> maptime = new HashMap<Long, String>();
		for (GPXEntry g : list_pre_my_gpx1) {
			maptime.put(g.getTime(), "");
		}

		// 
		MatchResult mr = mapMatching_my.doWork(list_pre_my_gpx1);

		List<PointOrList> result = new ArrayList<PointOrList>();

		int start = 0;
		int end = 0;
		GPXEntry startGPXEntry = null;
		GPXEntry endGPXEntry = null;
		boolean bool = true;
		AlternativeRoute altDijkstra = new AlternativeRoute(ghs, weighting,
				traversalMode);
		altDijkstra.setMaxShareFactor(0.7);
		altDijkstra.setMaxWeightFactor(3);
		altDijkstra.setMaxPaths(maxPath);
		/**
		 * 
		 * 
		 * 
		 * inputGPXEntries
		 */
		int count = 0;
		for (int i = 0; i < mr.getEdgeMatches().size(); i++) {

			int size = mr.getEdgeMatches().get(i).getGpxExtensions().size();
			// snapped
			if (size > 0 && bool) {
				count++;
				start = i;
				startGPXEntry = mr.getEdgeMatches().get(i).getGpxExtensions()
						.get(size - 1).getEntry();

				bool = false;
			} else if (!bool && size > 0) {
				end = i;
				Dijkstra dij = new Dijkstra(ghs, weighting, traversalMode);
				double d = dij.calcPath(
						mr.getEdgeMatches().get(start).getGpxExtensions()
								.get(0).getQueryResult().getClosestEdge()
								.getBaseNode(),
						mr.getEdgeMatches().get(end).getGpxExtensions().get(0)
								.getQueryResult().getClosestEdge()
								.getBaseNode()).getDistance();
				endGPXEntry = mr.getEdgeMatches().get(i).getGpxExtensions()
						.get(0).getEntry();

				double disOber = DistanceTwoPoint.Distance(startGPXEntry.lat,
						startGPXEntry.lon, endGPXEntry.lat, endGPXEntry.lon);
				if (d / 1000 / disOber > 1.20) {
					// 
					count++;
					List<AlternativeRoute.AlternativeInfo> pathInfos = null;
					pathInfos = altDijkstra.calcAlternatives(mr
							.getEdgeMatches().get(start).getGpxExtensions()
							.get(0).getQueryResult().getClosestEdge()
							.getBaseNode(), mr.getEdgeMatches().get(end)
							.getGpxExtensions().get(0).getQueryResult()
							.getClosestEdge().getBaseNode());

					List<Route> pathInfos1 = estimateAlt(pathInfos,
							startGPXEntry, endGPXEntry, maptime);
					GPXEntry g = new GPXEntry(0, 0, 0);
					PointOrList pl = new PointOrList(pathInfos1, g);
					result.add(pl);

				} else {
					List<Route> pathInfos = null;

					List<GPXEntry> templist = new ArrayList<GPXEntry>();

					for (int x = start; x < end + 1; x++) {
						PointList pl = mr.getEdgeMatches().get(x)
								.getEdgeState()
								.fetchWayGeometry(x == 0 ? 3 : 2);
						for (int k = 0; k < pl.getSize(); k++) {
							GPXEntry g = new GPXEntry(pl.getLatitude(k),
									pl.getLongitude(k), 0);
							templist.add(g);
						}
					}
					double dist = calDist(templist);
					int all = 0;
					double speed = dist
							* 1000000
							/ (Math.abs(startGPXEntry.getTime()
									- endGPXEntry.getTime()));
					double lon = mr.getEdgeMatches().get(start).getEdgeState()
							.fetchWayGeometry(start == 0 ? 3 : 2)
							.getLongitude(0);
					double lat = mr.getEdgeMatches().get(start).getEdgeState()
							.fetchWayGeometry(start == 0 ? 3 : 2)
							.getLatitude(0);
					GPXEntry gstart = new GPXEntry(lat, lon,
							startGPXEntry.getTime());
					String gPXEntry = getLast(result);

					for (int x = start; x < end + 1; x++) {
						PointList pl = mr.getEdgeMatches().get(x)
								.getEdgeState()
								.fetchWayGeometry(x == 0 ? 3 : 2);
						if (x == start
								&& !gPXEntry.equals(pl.getLat(0) + ","
										+ pl.getLon(0)) && start != 0) {
							for (int k = 1; k < pl.size(); k++) {
								double dis = calDist1(templist, all);

								long tplus = (long) (dis * 1000000 / speed);
								long t = gstart.getTime() + tplus;
								if (maptime.containsKey(t)) {
									PointOrList pll = new PointOrList(
											pathInfos, new GPXEntry(
													pl.getLat(k), pl.getLon(k),
													1, t));
									result.add(pll);

									all++;
								} else {
									PointOrList pll = new PointOrList(
											pathInfos, new GPXEntry(
													pl.getLat(k), pl.getLon(k),
													0, t));
									result.add(pll);

									all++;
								}

							}
						} else {
							for (int k = 0; k < pl.size(); k++) {
								double dis = calDist1(templist, all);

								long tplus = (long) (dis * 1000000 / speed);
								long t = gstart.getTime() + tplus;
								if (maptime.containsKey(t)) {
									PointOrList pll = new PointOrList(
											pathInfos, new GPXEntry(
													pl.getLat(k), pl.getLon(k),
													1, t));
									result.add(pll);

									all++;
								} else {
									PointOrList pll = new PointOrList(
											pathInfos, new GPXEntry(
													pl.getLat(k), pl.getLon(k),
													0, t));
									result.add(pll);

									all++;
								}

							}
						}

					}
				}
				// startGPXEntry = endGPXEntry;
				startGPXEntry = mr.getEdgeMatches().get(i).getGpxExtensions()
						.get(size - 1).getEntry();
				start = i;
			}

		}
		return result;

	}

	private static List<Point> readLine(String str) {
		List<Point> list = new ArrayList<Point>();
		if (!(str.indexOf("@") == -1)) {

			String[] seg1 = str.split("@", -1);
			if (seg1.length > 1) {
				for (int i = 0; i < seg1.length; i++) {
					String[] ss = seg1[i].split(",");
					double lat = Double.parseDouble(ss[1]);
					double lon = Double.parseDouble(ss[2]);
					Point p = new Point(ss[0], lon, lat);
					list.add(p);
				}
			}
		}

		return list;
	}

	private static double calDist1(List<GPXEntry> templist, int count) {
		double dist = 0;
		for (int i = 0; i < count; i++) {
			dist = dist
					+ DistanceTwoPoint.Distance(templist.get(i).lat,
							templist.get(i).lon, templist.get(i + 1).lat,
							templist.get(i + 1).lon);
		}

		return dist;
	}

	private static String getLast(List<PointOrList> result) {
		if (result.size() > 0) {
			int len = result.size() - 1;
			if (result.get(len).getList() != null) {
				List<GPXEntry> list = result.get(len).getList().get(0)
						.getList();
				GPXEntry g = list.get(list.size() - 1);
				return g.lat + "," + g.lon;
			} else {
				return result.get(len).getP().lat + ","
						+ result.get(len).getP().lon;
			}
		} else {
			return "";
		}

	}

	private static List<Route> estimateAlt(List<AlternativeInfo> pathInfos,
			GPXEntry startGPXEntry, GPXEntry endGPXEntry,
			Map<Long, String> maptime) {
		List<Route> pathInfos1 = new ArrayList<Route>();
		for (int i = 0; i < pathInfos.size(); i++) {
			List<GPXEntry> list = new ArrayList<GPXEntry>();
			double speed = pathInfos.get(i).getPath().getDistance()
					* 1000
					/ (Math.abs(startGPXEntry.getTime() - endGPXEntry.getTime()));

			PointList pl = pathInfos.get(i).getPath().calcPoints();
			for (int j = 0; j < pl.size(); j++) {
				double dis = 1000 * calDistEdgeIteratorState(pl, j);

				long t = startGPXEntry.getTime() + (long) (dis / speed) * 1000;
				if (maptime.containsKey(t)) {
					GPXEntry gx = new GPXEntry(pl.getLatitude(j), pl.getLon(j),
							1, t);
					list.add(gx);
				} else {
					GPXEntry gx = new GPXEntry(pl.getLatitude(j), pl.getLon(j),
							0, t);
					list.add(gx);
				}

			}

			Route route = new Route(list);
			pathInfos1.add(route);
		}

		return pathInfos1;
	}


	public List<GPXEntry> transToFusingTraj(MatchResult mr,
			List<GPXEntry> inputGPXEntries) throws ParseException {
		List<GPXEntry> inputGPXEntries2 = new ArrayList<GPXEntry>();
		long time = 0;
		int count = 0;
		for (int emIndex = 0; emIndex < mr.getEdgeMatches().size(); emIndex++) {
			EdgeMatch em = mr.getEdgeMatches().get(emIndex);

			if (em.getGpxExtensions().size() > 0) {
				for (int i = 0; i < em.getGpxExtensions().size(); i++) {
					GHPoint3D g = em.getGpxExtensions().get(i).getQueryResult()
							.getSnappedPoint();
					GPXEntry ge = new GPXEntry(g.lat, g.lon, 0, inputGPXEntries
							.get(count).getTime());
					inputGPXEntries2.add(ge);
					count++;
				}

			} else {
				PointList pl = em.getEdgeState().fetchWayGeometry(
						emIndex == 0 ? 3 : 2);

				for (int i = 0; i < pl.size(); i++) {

					inputGPXEntries2.add(new GPXEntry(pl.getLatitude(i), pl
							.getLongitude(i), 1, time));
				}
				// count++;
			}

		}
		return inputGPXEntries2;
	}


	private static List<GPXEntry> PointtoGPXEntry(List<Point> list)
			throws ParseException {
		List<GPXEntry> result = new ArrayList<GPXEntry>();
		DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		for (int i = 0; i < list.size(); i++) {
			double lat = list.get(i).getLat();
			double lon = list.get(i).getLon();
			Date d = df.parse(list.get(i).getTime());
			GPXEntry gx = new GPXEntry(lat, lon, 0, d.getTime());
			result.add(gx);
		}
		return result;
	}

	private static double calDistEdgeIteratorState(PointList pl, int j) {
		double dist = 0;
		for (int i = 0; i < j; i++) {
			dist = dist
					+ DistanceTwoPoint.Distance(pl.getLatitude(i),
							pl.getLongitude(i), pl.getLatitude(i + 1),
							pl.getLongitude(i + 1));
		}

		return dist;
	}

	private static double calDist(List<GPXEntry> templist) {
		double dist = 0;
		for (int i = 0; i < templist.size() - 1; i++) {
			dist = dist
					+ DistanceTwoPoint.Distance(templist.get(i).lat,
							templist.get(i).lon, templist.get(i + 1).lat,
							templist.get(i + 1).lon);
		}

		return dist;
	}

	private static Map<String, String> readEdgeDensity(String path)
			throws IOException {
		Map<String, String> map = new HashMap<String, String>();
		File file = new File(path);
		String encoding = "GBK";
		if (file.isFile() && file.exists()) { //  
			InputStreamReader read1 = new InputStreamReader(
					new FileInputStream(path), encoding);//  
			BufferedReader bufferedReader1 = new BufferedReader(read1);
			String lineTxt1 = null;

			while ((lineTxt1 = bufferedReader1.readLine()) != null) {
				String[] s = lineTxt1.split("\t");
				String[] seg = s[0].split(",");
				map.put(seg[2] + "," + seg[3], s[1]);
			}
		}

		return map;
	}

	public static void main(String[] args) throws Exception {
		Configuration conf = new Configuration();
		conf.set("mapred.textoutputformat.ignoreseparator", "true");
		conf.set("mapred.textoutputformat.separator", "");
		conf.set("mapreduce.job.queuename", args[6]);
		// conf.set("mapreduce.job.priority", "VERY_HIGH");
		conf.set("mapred.child.java.opts", "-Xmx" + args[7] + "m");
		conf.set("mapreduce.map.failures.maxpercent", "1");
		conf.set("mapred.task.timeout", "0");
		conf.set("gpsAccuracy", args[3]);
		conf.set("maxNode", args[4]);
		conf.set("maxpath", args[5]);
		Path vocabularyPath = new Path(args[2]);
		DistributedCache.addCacheFile(vocabularyPath.toUri(), conf);
		Job job = new Job(conf, "MM");
		// job.addCacheFile(new Path(args[2]).toUri());
		job.setJarByClass(MM_my_hadoop_multi.class);
		job.setMapperClass(MMMap.class);
		job.setNumReduceTasks(0);
		job.setOutputValueClass(Text.class);
		FileInputFormat.addInputPath(job, new Path(args[0]));
		FileOutputFormat.setOutputPath(job, new Path(args[1]));
		System.exit(job.waitForCompletion(true) ? 0 : 1);
	}

}
