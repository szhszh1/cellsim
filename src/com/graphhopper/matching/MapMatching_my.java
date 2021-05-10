package com.graphhopper.matching;

import com.graphhopper.GraphHopper;
import com.graphhopper.matching.util.HmmProbabilities;
import com.graphhopper.matching.util.TimeStep;
import com.graphhopper.routing.weighting.Weighting;
import com.bmw.hmm.SequenceState;
import com.bmw.hmm.ViterbiAlgorithm;
import com.graphhopper.routing.AlgorithmOptions;
import com.graphhopper.routing.Path;
import com.graphhopper.routing.QueryGraph;
import com.graphhopper.routing.RoutingAlgorithm;
import com.graphhopper.routing.RoutingAlgorithmFactory;
import com.graphhopper.routing.ch.CHAlgoFactoryDecorator;
import com.graphhopper.routing.ch.PrepareContractionHierarchies;
import com.graphhopper.routing.util.*;
import com.graphhopper.routing.weighting.FastestWeighting;
import com.graphhopper.storage.CHGraph;
import com.graphhopper.storage.Graph;
import com.graphhopper.storage.index.LocationIndexTree;
import com.graphhopper.storage.index.QueryResult;
import com.graphhopper.util.*;
import com.graphhopper.util.shapes.GHPoint;
import com.xjtu.util.DistanceTwoPoint;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class MapMatching_my {

	private final Graph routingGraph;
	private final LocationIndexMatch locationIndex;
	private double measurementErrorSigma = 50.0;
	private double transitionProbabilityBeta = 0.00959442;
	// private double transitionProbabilityBeta = 0.015;
	private final int nodeCount;
	private DistanceCalc distanceCalc = new DistancePlaneProjection();
	private final RoutingAlgorithmFactory algoFactory;
	private final AlgorithmOptions algoOptions;

	public MapMatching_my(GraphHopper hopper, AlgorithmOptions algoOptions) {
		this.locationIndex = new LocationIndexMatch(
				hopper.getGraphHopperStorage(),
				(LocationIndexTree) hopper.getLocationIndex());

		// create hints from algoOptions, so we can create the algorithm factory
		HintsMap hints = new HintsMap();
		for (Entry<String, String> entry : algoOptions.getHints().toMap()
				.entrySet()) {
			hints.put(entry.getKey(), entry.getValue());
		}

		// default is non-CH
		if (!hints.has(Parameters.CH.DISABLE)) {
			hints.put(Parameters.CH.DISABLE, true);
		}

		String vehicle = hints.getVehicle();
		if (vehicle.isEmpty()) {
			if (algoOptions.hasWeighting()) {
				vehicle = algoOptions.getWeighting().getFlagEncoder()
						.toString();
			} else {
				vehicle = hopper.getEncodingManager().fetchEdgeEncoders()
						.get(0).toString();
			}
			hints.setVehicle(vehicle);
		}

		if (!hopper.getEncodingManager().supports(vehicle)) {
			throw new IllegalArgumentException("Vehicle " + vehicle
					+ " unsupported. " + "Supported are: "
					+ hopper.getEncodingManager());
		}

		algoFactory = hopper.getAlgorithmFactory(hints);

		Weighting weighting = null;
		CHAlgoFactoryDecorator chFactoryDecorator = hopper
				.getCHFactoryDecorator();
		boolean forceFlexibleMode = hints.getBool(Parameters.CH.DISABLE, false);
		if (chFactoryDecorator.isEnabled() && !forceFlexibleMode) {
			if (!(algoFactory instanceof PrepareContractionHierarchies)) {
				throw new IllegalStateException(
						"Although CH was enabled a non-CH algorithm factory was returned "
								+ algoFactory);
			}

			weighting = ((PrepareContractionHierarchies) algoFactory)
					.getWeighting();
			this.routingGraph = hopper.getGraphHopperStorage().getGraph(
					CHGraph.class, weighting);
		} else {
			weighting = algoOptions.hasWeighting() ? algoOptions.getWeighting()
					: new FastestWeighting(hopper.getEncodingManager()
							.getEncoder(vehicle), algoOptions.getHints());
			this.routingGraph = hopper.getGraphHopperStorage();
		}

		this.algoOptions = AlgorithmOptions.start(algoOptions)
				.weighting(weighting).build();
		this.nodeCount = routingGraph.getNodes();
	}

	public void setDistanceCalc(DistanceCalc distanceCalc) {
		this.distanceCalc = distanceCalc;
	}

	/**
	 * Beta parameter of the exponential distribution for modeling transition
	 * probabilities.
	 */
	public void setTransitionProbabilityBeta(double transitionProbabilityBeta) {
		this.transitionProbabilityBeta = transitionProbabilityBeta;
	}

	/**
	 * Standard deviation of the normal distribution [m] used for modeling the
	 * GPS error.
	 */
	public void setMeasurementErrorSigma(double measurementErrorSigma) {
		this.measurementErrorSigma = measurementErrorSigma;
	}

	/**
	 * This method does the actual map matching.
	 * <p>
	 * 
	 * @param gpxList
	 *            the input list with GPX points which should match to edges of
	 *            the graph specified in the constructor
	 */
	public MatchResult doWork(List<GPXEntry> gpxList) {
		if (gpxList.size() < 2) {
			throw new IllegalArgumentException(
					"Too few coordinates in input file (" + gpxList.size()
							+ "). Correct format?");
		}
		final EdgeFilter edgeFilter = new DefaultEdgeFilter(algoOptions
				.getWeighting().getFlagEncoder());

		final List<QueryResult> allCandidates = new ArrayList<>();
		List<TimeStep<GPXExtension, GPXEntry, Path>> timeSteps = createTimeSteps(
				gpxList, edgeFilter, allCandidates);
//		System.out.println(gpxList.size()+"      "+timeSteps.size());
		if (allCandidates.size() < 2) {
			throw new IllegalArgumentException("Too few matching coordinates ("
					+ allCandidates.size() + "). Wrong region imported?");
		}
		if (timeSteps.size() < 2) {
			throw new IllegalStateException(
					"Coordinates produced too few time steps "
							+ timeSteps.size() + ", gpxList:" + gpxList.size());
		}

		final QueryGraph queryGraph = new QueryGraph(routingGraph)
				.setUseEdgeExplorerCache(true);

		queryGraph.lookup(allCandidates);
		for(int i=0;i<timeSteps.size();i++){
			System.out.print(i+"  ");
			System.out.println(timeSteps.get(i).candidates.size());
		}
		List<SequenceState<GPXExtension, GPXEntry, Path>> seq = computeViterbiSequence(
				timeSteps, gpxList, queryGraph);
//		System.out.println(seq.get(1).state.getQueryResult().getSnappedPoint().toString());
//		DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		// for(int i=0;i<seq.size();i++){
		// int a = seq.get(i).state.getQueryResult().getClosestNode();
		// System.out.println(df.format(new
		// Date(gpxList.get(i).getTime()))+"\t"+queryGraph.getNodeAccess().getLongitude(a)+"\t"+queryGraph.getNodeAccess().getLatitude(a));
		// }
		final EdgeExplorer explorer = queryGraph.createEdgeExplorer(edgeFilter);

		MatchResult matchResult = computeMatchResult(seq, gpxList,
				allCandidates, explorer);

		return matchResult;
	}

	/**
	 * Creates TimeSteps for the GPX entries but does not create emission or
	 * transition probabilities.
	 * 
	 * @param outAllCandidates
	 *            output parameter for all candidates, must be an empty list.
	 */
	private List<TimeStep<GPXExtension, GPXEntry, Path>> createTimeSteps(
			List<GPXEntry> gpxList, EdgeFilter edgeFilter,
			List<QueryResult> outAllCandidates) {
		int indexGPX = 0;
		TimeStep<GPXExtension, GPXEntry, Path> prevTimeStep = null;
		final List<TimeStep<GPXExtension, GPXEntry, Path>> timeSteps = new ArrayList<>();

		for (GPXEntry gpxEntry : gpxList) {
//			String ll = gpxEntry.lat+","+gpxEntry.lon;
//			double acc = 200;
//			if(map.containsKey(ll)){
//				acc = getSigma(ll,map);
//			}
			if (prevTimeStep == null
					|| distanceCalc.calcDist(prevTimeStep.observation.getLat(),
							prevTimeStep.observation.getLon(),
							gpxEntry.getLat(), gpxEntry.getLon()) > 200
					// always include last point
					|| indexGPX == gpxList.size() - 1
					){

				final List<QueryResult> queryResults = locationIndex
						.findNClosest(gpxEntry.lat, gpxEntry.lon, edgeFilter,
								measurementErrorSigma);
//				System.out.println(gpxEntry.lat+"   "+ gpxEntry.lon);
				outAllCandidates.addAll(queryResults);
				final List<GPXExtension> candidates = new ArrayList<>();
				for (QueryResult candidate : queryResults) {
					candidates.add(new GPXExtension(gpxEntry, candidate,
							indexGPX));
				}
				final TimeStep<GPXExtension, GPXEntry, Path> timeStep = new TimeStep<>(
						gpxEntry, candidates);
				timeSteps.add(timeStep);
				prevTimeStep = timeStep;
			}
			indexGPX++;
		}
		return timeSteps;
	}

	private List<SequenceState<GPXExtension, GPXEntry, Path>> computeViterbiSequence(
			List<TimeStep<GPXExtension, GPXEntry, Path>> timeSteps,
			List<GPXEntry> gpxList, final QueryGraph queryGraph) {
		final HmmProbabilities probabilities = new HmmProbabilities(
				measurementErrorSigma, transitionProbabilityBeta);
		final ViterbiAlgorithm<GPXExtension, GPXEntry, Path> viterbi = new ViterbiAlgorithm<>();

		int timeStepCounter = 0;
		TimeStep<GPXExtension, GPXEntry, Path> prevTimeStep = null;
		
		List<Double> dir = new ArrayList<Double>();
		for (int i=0;i<timeSteps.size()-1;i++) {
			double dirOber = calDir(timeSteps.get(i).observation, timeSteps.get(i+1).observation);
//			System.out.println(dirOber);
			if(dirOber==0){
				dirOber=dirOber+0.001;
			}
			dir.add(dirOber);
			
		}
		dir.add(dir.get(dir.size()-1));
		
		for (TimeStep<GPXExtension, GPXEntry, Path> timeStep : timeSteps) {
			computeEmissionProbabilities(timeStep, probabilities,dir.get(timeStepCounter));

			if (prevTimeStep == null) {
				viterbi.startWithInitialObservation(timeStep.observation,
						timeStep.candidates, timeStep.emissionLogProbabilities);
			} else {
				computeTransitionProbabilities(prevTimeStep, timeStep,
						probabilities, queryGraph);
				viterbi.nextStep(timeStep.observation, timeStep.candidates,
						timeStep.emissionLogProbabilities,
						timeStep.transitionLogProbabilities, timeStep.roadPaths);
			}
			if (viterbi.isBroken()) {
				String likelyReasonStr = "";
				if (prevTimeStep != null) {
					GPXEntry prevGPXE = prevTimeStep.observation;
					GPXEntry gpxe = timeStep.observation;
					double dist = distanceCalc.calcDist(prevGPXE.lat,
							prevGPXE.lon, gpxe.lat, gpxe.lon);
					
					if (dist > 2000) {
						likelyReasonStr = "Too long distance to previous measurement? "
								+ Math.round(dist) + "m, ";
					}
				}

				throw new RuntimeException(
						"Sequence is broken for submitted track at time step "
								+ timeStepCounter
								+ " ("
								+ gpxList.size()
								+ " points). "
								+ likelyReasonStr
								+ "observation:"
								+ timeStep.observation
								+ ", "
								+ timeStep.candidates.size()
								+ " candidates: "
								+ getSnappedCandidates(timeStep.candidates)
								+ ". If a match is expected consider increasing max_visited_nodes.");
			}

			timeStepCounter++;
			prevTimeStep = timeStep;
		}

		return viterbi.computeMostLikelySequence();
	}

	private void computeEmissionProbabilities(
			TimeStep<GPXExtension, GPXEntry, Path> timeStep,
			HmmProbabilities probabilities, Double dir) {
		FlagEncoder f = algoOptions.getWeighting().getFlagEncoder();
		
		for (GPXExtension candidate : timeStep.candidates) {
			String ll = candidate.getEntry().lat+","+candidate.getEntry().lon;
			double dcand = calDir1(candidate.getQueryResult().getClosestEdge().fetchWayGeometry(3));
			double jiajiao = Math.abs(dcand-dir);
			double edge_speed = f.getSpeed(candidate.getQueryResult().getClosestEdge().getFlags());
			//0-180    5-120
//			double wr = computeHybridW(jiajiao,edge_speed);
			
			double wspeed = computeSpeedWr(edge_speed);
			double wdir = computeDirWr(jiajiao);
			final double distance = candidate.getQueryResult()
					.getQueryDistance() *10*wspeed*wdir;
//			double sigma = getSigma(ll,map);
			
			timeStep.addEmissionLogProbability(candidate,
					probabilities.emissionLogProbability(distance));
		}
	}
	private double getSigma(String ll, Map<String, String> map) {
		if(map.containsKey(ll)){
			
			String str = map.get(ll);
			String[] seg = str.split("@");
			for(int i=0;i<seg.length;i++){
				int edge = Integer.parseInt(seg[i]);
				if(edge>=10){
					return i*50;
				}
			}
		}
		return 500;
	}

	private double computeDirWr(double jiajiao) {
		
		return jiajiao/180;
	}

	private double computeSpeedWr(double speed) {
		double wr = 1;
		if (speed < 31) {
			wr=1;
		} else{
			wr= 1- (0.09* ((speed / 10) - 1));
//			wr=  1-(0.08* ((speed / 10) - 1));
		}
		return wr;
	}
	private void computeTransitionProbabilities(
			TimeStep<GPXExtension, GPXEntry, Path> prevTimeStep,
			TimeStep<GPXExtension, GPXEntry, Path> timeStep,
			HmmProbabilities probabilities, QueryGraph queryGraph) {
		Map<String,Path> map = new HashMap<String,Path>();
		double minDistance=Double.MAX_VALUE;
		//实现最短距离的计算
		for (GPXExtension from : prevTimeStep.candidates) {
			for (GPXExtension to : timeStep.candidates) {
//				if(contains(from,prevTimeStep) && contains(to,timeStep)){
					RoutingAlgorithm algo1 = algoFactory.createAlgo(queryGraph,algoOptions);
					final Path path = algo1.calcPath(from.getQueryResult().getClosestNode(), to.getQueryResult().getClosestNode());
					if (path.isFound()) {
						map.put(from.toString()+","+to.toString(),path);
						if(minDistance > path.getDistance()){
							minDistance = path.getDistance();
						}
					}
//				}
//				else{
//					RoutingAlgorithm algo1 = algoFactory.createAlgo(queryGraph,algoOptions);
//					final Path path = algo1.calcPath(from.getQueryResult().getClosestNode(), to.getQueryResult().getClosestNode());
//					if (path.isFound()) {
//						if(minDistance>path.getDistance()*3){
//							minDistance=path.getDistance()*3;
//						}
//					}
//				}
			}
		}
		
//		System.out.println(minDistance);
		final double timeDiff = (timeStep.observation.getTime() - prevTimeStep.observation
				.getTime()) / 1000.0;
		for (GPXExtension from : prevTimeStep.candidates) {
			for (GPXExtension to : timeStep.candidates) {
				if (map.containsKey(from.toString()+","+to.toString())){
					Path path = map.get(from.toString()+","+to.toString());
					timeStep.addRoadPath(from, to, path);
//				RoutingAlgorithm algo = algoFactory.createAlgo(queryGraph,algoOptions);
//				final Path path = algo.calcPath(from.getQueryResult().getClosestNode(), to.getQueryResult().getClosestNode());
				if (path.getDistance()!=0) {
//					timeStep.addRoadPath(from, to, path);
//					int wei = getPathWeight(path);
//					if(from.getQueryResult().getClosestEdge().getName().equals(to.getQueryResult().getClosestEdge().getName())){
						
					final double transitionLogProbability = probabilities.transitionLogProbability(path.getDistance(),
										minDistance, timeDiff);
//						System.out.println(from.getQueryResult().getClosestEdge().getName()+"   "+to.getQueryResult().getClosestEdge().getName()+"   "+transitionLogProbability);
					timeStep.addTransitionLogProbability(from, to, transitionLogProbability);
				}
//					else{
//						final double transitionLogProbability = probabilities
//								.transitionLogProbability(path.getDistance(),
//										minDistance, timeDiff);
//
//						timeStep.addTransitionLogProbability(from, to,
//								transitionLogProbability);
//					}
					

				}
			}
		}
//		}
		
	}

	private int getPathWeight(Path path) {
		int count = 1;
		for(int i=0;i<path.calcEdges().size()-1;i++){
			double d1 = calDir1(path.calcEdges().get(i).fetchWayGeometry(3));
			double d2 = calDir1(path.calcEdges().get(i+1).fetchWayGeometry(3));
			if(Math.abs(d1-d2)>75){
				count++;
			}
		}
		return count;
	}

	private double computeWr1(double speed) {
		if(speed<31){
			return 100;
		}
		
		return 0.001;
	}

	private boolean contains(GPXExtension to,
			TimeStep<GPXExtension, GPXEntry, Path> timeStep) {
		double lat = timeStep.observation.lat;
		double lon = timeStep.observation.lon;
		PointList pl = to.getQueryResult().getClosestEdge().fetchWayGeometry(3);
		double startLat = pl.getLatitude(0);
		double startLon = pl.getLongitude(0);
		double endLat = pl.getLatitude(pl.size()-1);
		double endLon = pl.getLongitude(pl.size()-1);
		if((startLat <= lat && endLat > lat) || (startLon <=lon && endLon >lon) || (startLat >= lat && endLat < lat) || (startLon >=lon && endLon <lon)){
			return true;
		}
		return false;
	}

	private double calWp(GPXExtension from, GPXExtension to) {
		if (from.getQueryResult().getClosestEdge().getEdge() == to
				.getQueryResult().getClosestEdge().getEdge()) {
			return 0.75;
		} else {
			return 0.25;
		}

	}
	private static double calDir(GPXEntry a, GPXEntry b) {
		double angle1 = 0;
		double lat1 = a.getLat();
		double lon1 = a.getLon();

		double lat2 = b.getLat();
		double lon2 = b.getLon();
		double lat3 = lat1;
		double lon3 = lon2;
		double dis12 = DistanceTwoPoint.Distance(lat1, lon1, lat2, lon2);
		double dis13 = DistanceTwoPoint.Distance(lat1, lon1, lat3, lon3);

		if (Math.abs(dis12 - dis13) < 0.001) {
			dis12 = dis13;
		}

		if (lon2 < lon1 && lat2 > lat1) {
			angle1 = 180 - Math.acos(dis13 / dis12) / Math.PI * 180;
		}

		if (lon2 < lon1 && lat2 < lat1) {
			angle1 = 180 + Math.acos(dis13 / dis12) / Math.PI * 180;
		}

		if (lon2 > lon1 && lat2 < lat1) {
			angle1 = 360 - Math.acos(dis13 / dis12) / Math.PI * 180;
		}

		if (lon2 > lon1 && lat2 > lat1) {
			angle1 = Math.acos(dis13 / dis12) / Math.PI * 180;
		}

		if (lon2 == lon1 && lat2 > lat1) {
			angle1 = 90;
		}

		if (lon2 == lon1 && lat2 < lat1) {
			angle1 = 270;
		}

		if (lon2 > lon1 && lat2 == lat1) {
			angle1 = 0;
		}

		if (lon2 < lon1 && lat2 == lat1) {
			angle1 = 180;
		}

		// return 0;
		if (angle1 >= 180) {
			angle1 = angle1 - 360;
		}
		return angle1;

	}
	private double calDir1(PointList plist) {
		GPXEntry a = new GPXEntry(plist.getLatitude(0), plist.getLon(0), 0);
		GPXEntry b = new GPXEntry(plist.getLatitude(plist.size() - 1),
				plist.getLon(plist.size() - 1), 0);
		double d = calDir(a, b);
		return d;
	}

	private MatchResult computeMatchResult(
			List<SequenceState<GPXExtension, GPXEntry, Path>> seq,
			List<GPXEntry> gpxList, List<QueryResult> allCandidates,
			EdgeExplorer explorer) {
		final Map<String, EdgeIteratorState> virtualEdgesMap = new HashMap<>();
		for (QueryResult candidate : allCandidates) {
			fillVirtualEdges(virtualEdgesMap, explorer, candidate);
		}

		MatchResult matchResult = computeMatchedEdges(seq, virtualEdgesMap);

		// for (int emIndex = 0; emIndex < matchResult.getEdgeMatches().size();
		// emIndex++) {
		// EdgeMatch em = matchResult.getEdgeMatches().get(emIndex);
		// PointList pl = em.getEdgeState().fetchWayGeometry(emIndex == 0 ? 3 :
		// 2);
		// System.out.println(pl.toString());
		// }
		// 
		// computeGpxStats(gpxList, matchResult);

		return matchResult;
	}

	private MatchResult computeMatchedEdges(
			List<SequenceState<GPXExtension, GPXEntry, Path>> seq,
			Map<String, EdgeIteratorState> virtualEdgesMap) {
		List<EdgeMatch> edgeMatches = new ArrayList<>();
		double distance = 0.0;
		long time = 0;
		EdgeIteratorState currentEdge = null;
		List<GPXExtension> gpxExtensions = new ArrayList<>();
		GPXExtension queryResult = seq.get(0).state;
		gpxExtensions.add(queryResult);

		for (int j = 1; j < seq.size(); j++) {
			queryResult = seq.get(j).state;
			Path path = seq.get(j).transitionDescriptor;
			distance += path.getDistance();

			time += path.getTime();
			for (EdgeIteratorState edgeIteratorState : path.calcEdges()) {
				EdgeIteratorState directedRealEdge = resolveToRealEdge(
						virtualEdgesMap, edgeIteratorState);
				if (directedRealEdge == null) {
					throw new RuntimeException("Did not find real edge for "
							+ edgeIteratorState.getEdge());
				}
				if (currentEdge == null
						|| !equalEdges(directedRealEdge, currentEdge)) {
					if (currentEdge != null) {

						EdgeMatch edgeMatch = new EdgeMatch(currentEdge,
								gpxExtensions);
						edgeMatches.add(edgeMatch);
						gpxExtensions = new ArrayList<>();
					}
					currentEdge = directedRealEdge;
				}
			}
			// System.out.println(queryResult.queryResult.getSnappedPoint());

			gpxExtensions.add(queryResult);
			// System.out.println(gpxExtensions);
		}
		EdgeMatch lastEdgeMatch = edgeMatches.get(edgeMatches.size() - 1);

		if (!gpxExtensions.isEmpty()
				&& !equalEdges(currentEdge, lastEdgeMatch.getEdgeState())) {
			edgeMatches.add(new EdgeMatch(currentEdge, gpxExtensions));
		} else {
			lastEdgeMatch.getGpxExtensions().addAll(gpxExtensions);
		}
		MatchResult matchResult = new MatchResult(edgeMatches);

		matchResult.setMatchMillis(time);
		matchResult.setMatchLength(distance);
		// for(int i=0;i<matchResult.getEdgeMatches().size();i++){
		// System.out.println(matchResult.getEdgeMatches().get(i).getEdgeState());
		// }
		// for (final EdgeMatch match : edgeMatches) {
		// final EdgeIteratorState s = match.getEdgeState();
		// final int eId = s.getEdge();
		// System.out.print(eId);
		// System.out.print("\t");
		// System.out.println(s.fetchWayGeometry(1));
		// }

		return matchResult;
	}

	/**
	 * Calculate GPX stats to determine quality of matching.
	 */
	private void computeGpxStats(List<GPXEntry> gpxList, MatchResult matchResult) {
		double gpxLength = 0;
		GPXEntry prevEntry = gpxList.get(0);
		for (int i = 1; i < gpxList.size(); i++) {
			GPXEntry entry = gpxList.get(i);
			gpxLength += distanceCalc.calcDist(prevEntry.lat, prevEntry.lon,
					entry.lat, entry.lon);
			prevEntry = entry;
		}

		long gpxMillis = gpxList.get(gpxList.size() - 1).getTime()
				- gpxList.get(0).getTime();

		matchResult.setGPXEntriesMillis(gpxMillis);
		matchResult.setGPXEntriesLength(gpxLength);
	}

	private boolean equalEdges(EdgeIteratorState edge1, EdgeIteratorState edge2) {
		return edge1.getEdge() == edge2.getEdge()
				&& edge1.getBaseNode() == edge2.getBaseNode()
				&& edge1.getAdjNode() == edge2.getAdjNode();
	}

	private EdgeIteratorState resolveToRealEdge(
			Map<String, EdgeIteratorState> virtualEdgesMap,
			EdgeIteratorState edgeIteratorState) {
		if (isVirtualNode(edgeIteratorState.getBaseNode())
				|| isVirtualNode(edgeIteratorState.getAdjNode())) {
			return virtualEdgesMap.get(virtualEdgesMapKey(edgeIteratorState));
		} else {
			return edgeIteratorState;
		}
	}

	private boolean isVirtualNode(int node) {
		return node >= nodeCount;
	}

	/**
	 * Fills the minFactorMap with weights for the virtual edges.
	 */
	private void fillVirtualEdges(
			Map<String, EdgeIteratorState> virtualEdgesMap,
			EdgeExplorer explorer, QueryResult qr) {
		if (isVirtualNode(qr.getClosestNode())) {
			EdgeIterator iter = explorer.setBaseNode(qr.getClosestNode());
			while (iter.next()) {
				int node = traverseToClosestRealAdj(explorer, iter);
				if (node == qr.getClosestEdge().getAdjNode()) {
					virtualEdgesMap.put(virtualEdgesMapKey(iter), qr
							.getClosestEdge().detach(false));
					virtualEdgesMap.put(reverseVirtualEdgesMapKey(iter), qr
							.getClosestEdge().detach(true));
				} else if (node == qr.getClosestEdge().getBaseNode()) {
					virtualEdgesMap.put(virtualEdgesMapKey(iter), qr
							.getClosestEdge().detach(true));
					virtualEdgesMap.put(reverseVirtualEdgesMapKey(iter), qr
							.getClosestEdge().detach(false));
				} else {
					throw new RuntimeException();
				}
			}
		}
	}

	private String virtualEdgesMapKey(EdgeIteratorState iter) {
		return iter.getBaseNode() + "-" + iter.getEdge() + "-"
				+ iter.getAdjNode();
	}

	private String reverseVirtualEdgesMapKey(EdgeIteratorState iter) {
		return iter.getAdjNode() + "-" + iter.getEdge() + "-"
				+ iter.getBaseNode();
	}

	private int traverseToClosestRealAdj(EdgeExplorer explorer,
			EdgeIteratorState edge) {
		if (!isVirtualNode(edge.getAdjNode())) {
			return edge.getAdjNode();
		}

		EdgeIterator iter = explorer.setBaseNode(edge.getAdjNode());
		while (iter.next()) {
			if (iter.getAdjNode() != edge.getBaseNode()) {
				return traverseToClosestRealAdj(explorer, iter);
			}
		}
		throw new IllegalStateException("Cannot find adjacent edge " + edge);
	}

	private String getSnappedCandidates(Collection<GPXExtension> candidates) {
		String str = "";
		for (GPXExtension gpxe : candidates) {
			if (!str.isEmpty()) {
				str += ", ";
			}
			str += "distance: " + gpxe.queryResult.getQueryDistance() + " to "
					+ gpxe.queryResult.getSnappedPoint();
		}
		return "[" + str + "]";
	}

	private void printMinDistances(
			List<TimeStep<GPXExtension, GPXEntry, Path>> timeSteps) {
		TimeStep<GPXExtension, GPXEntry, Path> prevStep = null;
		int index = 0;
		for (TimeStep<GPXExtension, GPXEntry, Path> ts : timeSteps) {
			if (prevStep != null) {
				double dist = distanceCalc.calcDist(prevStep.observation.lat,
						prevStep.observation.lon, ts.observation.lat,
						ts.observation.lon);
				double minCand = Double.POSITIVE_INFINITY;
				for (GPXExtension prevGPXE : prevStep.candidates) {
					for (GPXExtension gpxe : ts.candidates) {
						GHPoint psp = prevGPXE.queryResult.getSnappedPoint();
						GHPoint sp = gpxe.queryResult.getSnappedPoint();
						double tmpDist = distanceCalc.calcDist(psp.lat,
								psp.lon, sp.lat, sp.lon);
						if (tmpDist < minCand) {
							minCand = tmpDist;
						}
					}
				}
				System.out.println(index + ": " + Math.round(dist)
						+ "m, minimum candidate: " + Math.round(minCand) + "m");
				index++;
			}

			prevStep = ts;
		}
	}

	// TODO: Make setFromNode and processEdge public in Path and then remove
	// this.
	private static class MyPath extends Path {

		public MyPath(Graph graph, Weighting weighting) {
			super(graph, weighting);
		}

		@Override
		public Path setFromNode(int from) {
			return super.setFromNode(from);
		}

		@Override
		public void processEdge(int edgeId, int adjNode, int prevEdgeId) {
			super.processEdge(edgeId, adjNode, prevEdgeId);
		}
	}

	public Path calcPath(MatchResult mr) {
		MyPath p = new MyPath(routingGraph, algoOptions.getWeighting());
		if (!mr.getEdgeMatches().isEmpty()) {
			int prevEdge = EdgeIterator.NO_EDGE;
			p.setFromNode(mr.getEdgeMatches().get(0).getEdgeState()
					.getBaseNode());
			for (EdgeMatch em : mr.getEdgeMatches()) {
				p.processEdge(em.getEdgeState().getEdge(), em.getEdgeState()
						.getAdjNode(), prevEdge);
				prevEdge = em.getEdgeState().getEdge();
			}

			// TODO p.setWeight(weight);
			p.setFound(true);

			return p;
		} else {
			return p;
		}
	}
}
