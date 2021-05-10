package com.xjtu.mm;

import java.io.IOException;
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

/**
 * 
 * @param input *** 
 * 		  output ***
 * 		  graph graph-cache 
 * 		  gpsAccuracy 200 
 * 		  column 3 
 * 		  maxNode 20000
 * 
 */
public class MM_my_hadoop_single {
	public static class MMMap extends Mapper<LongWritable, Text, Text, Text> {
		GraphHopper hopper;
		DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

		protected void setup(Context context) throws IOException,
				InterruptedException {
			Path[] localCacheFiles = DistributedCache.getLocalCacheFiles(context.getConfiguration());
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

				// my mapmatching
				MapMatching_my mapMatching_my = new MapMatching_my(
						hopper, opts);
				mapMatching_my.setTransitionProbabilityBeta(0.0959442);
				mapMatching_my.setMeasurementErrorSigma(gpsAccuracy);

				GraphHopperStorage ghs = hopper.getGraphHopperStorage();
				Weighting weighting = opts.getWeighting();
				TraversalMode traversalMode = opts.getTraversalMode();
				String[] sg = value.toString().split("\t");
				StringBuilder sb = new StringBuilder(sg[0]);
				sb.append("\t");
				/**
				 * Part 3  
				 */
				DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
				List<GPXEntry> result_cellular_my = mm_cellular_my(sg[1],
						mapMatching_my, ghs, weighting, traversalMode);
				for (int i = 0; i < result_cellular_my.size() - 1; i++) {
					sb.append(df.format(new Date(result_cellular_my.get(i).getTime())))
							.append(",")
							.append(String.valueOf(result_cellular_my.get(i).lat))
							.append(",")
							.append(String.valueOf(result_cellular_my.get(i).lon))
							.append(",");
					if (result_cellular_my.get(i).ele == 1) {
						sb.append("0");
					}
					if (result_cellular_my.get(i).ele == 0) {
						sb.append("1");
					}
					sb.append("&");
				}
				int len = result_cellular_my.size() - 1;
				sb.append(
						df.format(new Date(result_cellular_my.get(len)
								.getTime())))
						.append(",")
						.append(String.valueOf(result_cellular_my.get(len).lat))
						.append(",")
						.append(String.valueOf(result_cellular_my.get(len).lon))
						.append(",");
				if (result_cellular_my.get(len).ele == 1) {
					sb.append("0");
				}
				if (result_cellular_my.get(len).ele == 0) {
					sb.append("1");
				}

				word.set(sb.toString());
				val.set("");
				context.write(word, val);

			} catch (Exception e) {
			}
		}

	}

	private static List<GPXEntry> mm_cellular_my(String cellPath,
			MapMatching_my mapMatching_my, GraphHopperStorage ghs,
			Weighting weighting, TraversalMode traversalMode)
			throws ParseException, NumberFormatException, IOException {

		List<Point> list_my = readLine(cellPath);
		List<GPXEntry> list_pre_my_gpx1 = null;
		try {
			list_pre_my_gpx1 = PointtoGPXEntry(list_my);
		} catch (Exception e) {
		}
		Map<Long, String> maptime = new HashMap<Long, String>();
		List<GPXEntry> fusing = new ArrayList<GPXEntry>();
		for (GPXEntry g : list_pre_my_gpx1) {
			maptime.put(g.getTime(), "");
		}

		// 
		MatchResult mr = mapMatching_my.doWork(list_pre_my_gpx1);

		List<GPXEntry> fusingtraj = transToFusingTraj(mr, list_pre_my_gpx1);
		fusing = estimateTime(fusingtraj);
		return fusing;
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

	public static List<GPXEntry> transToFusingTraj(MatchResult mr,
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
			}

		}
		return inputGPXEntries2;
	}

	private static List<GPXEntry> estimateTime(List<GPXEntry> fusingtraj)
			throws ParseException {
		List<GPXEntry> result = new ArrayList<GPXEntry>();
		for (int i = 0; i < fusingtraj.size() - 1; i++) {
			double ele1 = fusingtraj.get(i).ele;

			if (ele1 == 0) {
				for (int j = i + 1; j < fusingtraj.size(); j++) {
					double ele2 = fusingtraj.get(j).ele;
					if (ele2 == 0) {

						List<GPXEntry> res = calTime(i, j, fusingtraj);
						result.add(fusingtraj.get(i));
						result.addAll(res);

						i = j - 1;
						break;
					}
				}

			}

		}
		result.add(fusingtraj.get(fusingtraj.size() - 1));
		return result;
	}

	private static List<GPXEntry> calTime(int i, int j,
			List<GPXEntry> fusingtraj) throws ParseException {
		List<GPXEntry> result = new ArrayList<GPXEntry>();
		double dis = 0;
		for (int k = i; k < j; k++) {
			dis += 1000 * DistanceTwoPoint.Distance(fusingtraj.get(k).lat,
					fusingtraj.get(k).lon, fusingtraj.get(k + 1).lat,
					fusingtraj.get(k + 1).lon);

		}
		long minus = fusingtraj.get(j).getTime() - fusingtraj.get(i).getTime();
		double speed = 0;
		if (minus != 0) {
			speed = dis / (minus / 1000);
		} else {
			speed = 0;
		}
		for (int k = i + 1; k < j; k++) {
			GPXEntry g = determineLoc(speed, i, k, fusingtraj);
			result.add(g);
		}

		return result;
	}

	private static GPXEntry determineLoc(double speed, int start, int end,
			List<GPXEntry> fusingtraj) {
		double disab = 0;
		for (int i = start; i <= end - 1; i++) {
			disab += 1000 * DistanceTwoPoint.Distance(fusingtraj.get(i)
					.getLat(), fusingtraj.get(i).getLon(), fusingtraj
					.get(i + 1).getLat(), fusingtraj.get(i + 1).getLon());
		}

		long time = (long) (disab / speed);
		GPXEntry ge = new GPXEntry(fusingtraj.get(end).lat,
				fusingtraj.get(end).lon, 1, fusingtraj.get(start).getTime()
						+ time * 1000);

		return ge;
	}

	public static void main(String[] args) throws Exception {
		Configuration conf = new Configuration();
		conf.set("mapred.textoutputformat.ignoreseparator", "true");
		conf.set("mapred.textoutputformat.separator", "");
		conf.set("mapreduce.job.queuename", args[5]);
		conf.set("mapreduce.job.priority", "VERY_HIGH");
		conf.set("mapred.child.java.opts", "-Xmx" + args[6] + "m");
		conf.set("mapreduce.map.failures.maxpercent", "1");
		conf.set("mapred.task.timeout", "0");
		conf.set("gpsAccuracy", args[3]);
		conf.set("maxNode", args[4]);
		Path vocabularyPath = new Path(args[2]);
		DistributedCache.addCacheFile(vocabularyPath.toUri(), conf);
		Job job = new Job(conf, "MM");
		job.setJarByClass(MM_my_hadoop_single.class);
		job.setMapperClass(MMMap.class);
		job.setNumReduceTasks(0);
		job.setOutputValueClass(Text.class);
		FileInputFormat.addInputPath(job, new Path(args[0]));
		FileOutputFormat.setOutputPath(job, new Path(args[1]));
		System.exit(job.waitForCompletion(true) ? 0 : 1);
	}

}
