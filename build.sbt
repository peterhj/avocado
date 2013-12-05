import AssemblyKeys._

assemblySettings

name := "Avocado"

version := "0.0.1"

organization := "edu.berkeley.cs.amplab"

scalaVersion := "2.9.3"

classpathTypes ++= Set ("orbit")

net.virtualvoid.sbt.graph.Plugin.graphSettings

libraryDependencies ++= Seq(
  "org.apache.spark" %% "spark-core" % "0.8.0-incubating",
  "org.scalatest" %% "scalatest" % "1.9.1" % "test",
  "edu.berkeley.cs.amplab.adam" % "adam-format" % "0.5.0-SNAPSHOT",
  "edu.berkeley.cs.amplab.adam" % "adam-commands" % "0.5.0-SNAPSHOT",
  "args4j" % "args4j" % "2.0.23",
  "commons-lang" % "commons-lang" % "2.6",
  "variant" % "variant" % "1.93",
  "tribble" % "tribble" % "1.93"
)

libraryDependencies ++= Seq(
  "org.eclipse.jetty.orbit" % "javax.servlet" % "2.5.0.v201103041518" artifacts (Artifact("javax.servlet", "jar", "jar")
  )
)

resolvers ++= Seq(
  "Typesafe" at "http://repo.typesafe.com/typesafe/releases",
  "Scala Tools Snapshots" at "http://scala-tools.org/repo-snapshots/",
  "Spray" at "http://repo.spray.cc",
  "Local" at "file:///Users/fnothaft/.m2/repository",
  "massie-maven" at "http://www.cs.berkeley.edu/~massie/maven/",
  "apache" at "https://repository.apache.org/content/repositories/releases",
  "hadoop-bam" at "http://hadoop-bam.sourceforge.net/maven/"
)

mergeStrategy in assembly := {
 case m if m.toLowerCase.endsWith("manifest.mf") => MergeStrategy.discard
 case m if m.toLowerCase.matches("meta-inf.*\\.sf$") => MergeStrategy.discard
 case "META-INF/services/org.apache.hadoop.fs.FileSystem" =>
MergeStrategy.concat
 case "reference.conf" => MergeStrategy.concat
 case "log4j.properties" => MergeStrategy.concat
 case _ => MergeStrategy.first
}

assemblyCacheOutput in assembly := false