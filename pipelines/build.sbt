assemblyJarName in assembly := "fgstr-pipelines-" + version.value + ".jar"
mainClass       in assembly := Some("dagr.core.cmdline.DagrCoreMain")

