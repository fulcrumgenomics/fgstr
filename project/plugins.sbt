resolvers += Resolver.url("fix-sbt-plugin-releases", url("http://dl.bintray.com/sbt/sbt-plugin-releases"))(Resolver.ivyStylePatterns)

addSbtPlugin("com.typesafe.sbt"  %  "sbt-git"        %  "0.8.5")
addSbtPlugin("com.eed3si9n"      %  "sbt-assembly"   %  "0.14.3")
addSbtPlugin("org.scoverage"     %  "sbt-scoverage"  %  "1.5.0")
addSbtPlugin("com.eed3si9n"      %  "sbt-unidoc"     %  "0.3.3")
