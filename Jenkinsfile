pipeline {
  parameters{string(defaultValue: '', description: 'api-token', name: 'API_TOKEN')
  	     string(defaultValue: '', description: 'Token for readthedocs', name: 'RTD_TOKEN')
	     string(defaultValue: '', description: 'buildable phid', name: 'BUILD_TARGET_PHID')
	     string(defaultValue: '', description: 'Commit id', name: 'COMMIT_ID')
	     string(defaultValue: '', description: 'Diff id', name: 'DIFF_ID')
	     string(defaultValue: 'PHID-PROJ-qr4btg2lycr6fpvkkmb5', description: 'ID of the project', name: 'PROJECT_ID')}
 
  environment{
    PHABRICATOR_HOST = 'https://c4science.ch/api/'
    PYTHONPATH = sh returnStdout: true, script: 'echo ${WORKSPACE}/tests/ci/script/'
  }

  agent {
    dockerfile { additionalBuildArgs '--tag cracklet-environment' }
  }
  options{
	skipDefaultCheckout(true)
  }

  stages {
    stage('SCM Checkout') {
      steps {
      	    // Clean before build
	    cleanWs()
	    // Checkout from SCM
      	    checkout scm: [
		 $class:'GitSCM',
		 branches:scm.branches,
		 extensions:[[
			$class:'SubmoduleOption',
		 	recursiveSubmodules:true,
		 ]],
		 userRemoteConfigs:scm.userRemoteConfigs]
		 }
	}
	}
post{
	always{
	sh"docker system prune -f"
	}
}
}